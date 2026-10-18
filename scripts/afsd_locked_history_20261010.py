"""Process-lifetime Windows read locks for immutable simulation evidence.

Cache authority is a live, read-only Windows file handle, not timestamps or a
receipt left by an earlier process. Files are opened with FILE_SHARE_READ only;
ancestor directory handles allow reads/writes but deny deletion/rename. SHA256
is computed after taking the locks, and every reuse checks the live handle and
the file identity obtained by reopening the original path. New output files can
still be created in locked directories. All locks disappear on process exit.

This protects normal filesystem writers/replacements, not privileged raw-disk
access or writable memory mappings established before the session. Only closed,
completed evidence files belong here; live logs and active outputs must not be
sealed. Final auditing rereads every byte, including files whose digest was
requested with sha() before an expected digest was supplied.
"""
from __future__ import annotations

import ctypes
from ctypes import wintypes
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import re
import threading


class HistoryLockError(RuntimeError):
    """Evidence cannot be kept immutable or does not match its expectation."""


class _FileInformation(ctypes.Structure):
    _fields_ = [
        ('attributes', wintypes.DWORD),
        ('creation', wintypes.FILETIME),
        ('access', wintypes.FILETIME),
        ('write', wintypes.FILETIME),
        ('volume', wintypes.DWORD),
        ('size_high', wintypes.DWORD),
        ('size_low', wintypes.DWORD),
        ('links', wintypes.DWORD),
        ('index_high', wintypes.DWORD),
        ('index_low', wintypes.DWORD),
    ]


class _Windows:
    GENERIC_READ = 0x80000000
    FILE_READ_ATTRIBUTES = 0x80
    FILE_SHARE_READ = 1
    FILE_SHARE_WRITE = 2
    OPEN_EXISTING = 3
    FILE_FLAG_BACKUP_SEMANTICS = 0x02000000
    FILE_FLAG_OPEN_REPARSE_POINT = 0x00200000
    FILE_FLAG_SEQUENTIAL_SCAN = 0x08000000
    FILE_ATTRIBUTE_DIRECTORY = 0x10
    FILE_ATTRIBUTE_REPARSE_POINT = 0x400

    def __init__(self):
        if os.name != 'nt':
            raise HistoryLockError('locked history requires Windows file sharing semantics')
        self.kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        self.create = self.kernel.CreateFileW
        self.create.argtypes = [wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD,
                                wintypes.LPVOID, wintypes.DWORD, wintypes.DWORD,
                                wintypes.HANDLE]
        self.create.restype = wintypes.HANDLE
        self.info = self.kernel.GetFileInformationByHandle
        self.info.argtypes = [wintypes.HANDLE, ctypes.POINTER(_FileInformation)]
        self.info.restype = wintypes.BOOL
        self.close = self.kernel.CloseHandle
        self.close.argtypes = [wintypes.HANDLE]
        self.close.restype = wintypes.BOOL
        self.seek = self.kernel.SetFilePointerEx
        self.seek.argtypes = [wintypes.HANDLE, ctypes.c_longlong,
                              ctypes.POINTER(ctypes.c_longlong), wintypes.DWORD]
        self.seek.restype = wintypes.BOOL
        self.read = self.kernel.ReadFile
        self.read.argtypes = [wintypes.HANDLE, wintypes.LPVOID, wintypes.DWORD,
                              ctypes.POINTER(wintypes.DWORD), wintypes.LPVOID]
        self.read.restype = wintypes.BOOL
        self.invalid = ctypes.c_void_p(-1).value

    @staticmethod
    def native_path(path):
        value = str(path)
        if value.startswith('\\\\'):
            return '\\\\?\\UNC\\' + value[2:]
        return '\\\\?\\' + value

    def open(self, path, directory=False):
        access = self.FILE_READ_ATTRIBUTES if directory else self.GENERIC_READ
        share = self.FILE_SHARE_READ | (self.FILE_SHARE_WRITE if directory else 0)
        flags = self.FILE_FLAG_OPEN_REPARSE_POINT
        flags |= self.FILE_FLAG_BACKUP_SEMANTICS if directory else self.FILE_FLAG_SEQUENTIAL_SCAN
        handle = self.create(self.native_path(path), access, share, None,
                             self.OPEN_EXISTING, flags, None)
        if handle == self.invalid:
            raise ctypes.WinError(ctypes.get_last_error())
        return handle

    def identity(self, handle, directory=False):
        row = _FileInformation()
        if not self.info(handle, ctypes.byref(row)):
            raise ctypes.WinError(ctypes.get_last_error())
        if row.attributes & self.FILE_ATTRIBUTE_REPARSE_POINT:
            raise HistoryLockError('reparse points/symlinks are not immutable evidence paths')
        if bool(row.attributes & self.FILE_ATTRIBUTE_DIRECTORY) != bool(directory):
            raise HistoryLockError('evidence path has an unexpected file type')
        # Size is useful as an additional assertion; dates are never authority.
        return (int(row.volume), (int(row.index_high) << 32) | int(row.index_low),
                (int(row.size_high) << 32) | int(row.size_low))


@dataclass
class _HeldFile:
    path: Path
    handle: int
    identity: tuple
    digest: str


class LockedHistory:
    """A fail-closed, nonpersistent SHA256 cache backed by live OS locks.

    seal(mapping) and verify(mapping) take root-relative or absolute paths and
    SHA256 hex digests. Both add newly encountered files and validate cache hits.
    sha(path) supplies the same lock-protected digest before an expectation is
    known. final_audit() rereads all sealed bytes and retains the locks; close()
    releases them. A failed operation invalidates the entire session.
    """

    def __init__(self, root):
        self.root = Path(os.path.abspath(os.fspath(root)))
        if not self.root.is_dir():
            raise HistoryLockError('history root must be an existing directory')
        self._win = _Windows()
        self._mutex = threading.RLock()
        self._files = {}
        self._directories = {}
        self._expected = {}
        self._closed = False
        self._failed = False
        self._read_buffer = ctypes.create_string_buffer(1024 * 1024)
        self._counts = dict(requests=0, hash_reads=0, hash_bytes=0, cache_hits=0,
                            identity_checks=0, unique_files=0, final_audits=0)

    def _ensure_open(self):
        if self._closed or self._failed:
            raise HistoryLockError('locked history session is closed or invalidated')

    def _path(self, value):
        raw = os.fspath(value)
        if not isinstance(raw, str) or not raw or '\x00' in raw:
            raise HistoryLockError('a nonempty normal path is required')
        # Reject aliases before Path/abspath could silently normalize them.
        if raw.startswith(('\\\\?\\', '\\\\.\\')):
            raise HistoryLockError('device paths are forbidden')
        drive, tail = os.path.splitdrive(raw)
        if drive and not os.path.isabs(raw):
            raise HistoryLockError('drive-relative paths are forbidden')
        for part in re.split(r'[\\/]', tail):
            if part in ('.', '..') or (part and part[-1] in (' ', '.')):
                raise HistoryLockError('noncanonical or escaping evidence path')
            if any(character in part for character in ':*?<>|"'):
                raise HistoryLockError('alternate streams and special paths are forbidden')
        path = Path(os.path.abspath(raw if os.path.isabs(raw) else self.root / raw))
        root_key, key = self._key(self.root), self._key(path)
        try:
            if os.path.commonpath([root_key, key]) != root_key or key == root_key:
                raise HistoryLockError('dependency outside project root or not a file')
        except ValueError as exc:
            raise HistoryLockError('dependency outside project root') from exc
        return path

    @staticmethod
    def _key(path):
        return os.path.normcase(str(path))

    def _lock_directories(self, path):
        # Lock all ancestors, including above root: otherwise renaming an outer
        # folder could change the meaning of an otherwise locked relative path.
        for directory in reversed(path.parents):
            key = self._key(directory)
            if key in self._directories:
                continue
            handle = self._win.open(directory, directory=True)
            try:
                identity = self._win.identity(handle, directory=True)
            except BaseException:
                self._win.close(handle)
                raise
            self._directories[key] = (handle, identity)

    def _check_identity(self, row):
        held = self._win.identity(row.handle)
        if held != row.identity:
            raise HistoryLockError('held evidence identity changed: ' + str(row.path))
        # Reopen the named file, rather than trusting path metadata or mtime.
        handle = self._win.open(row.path)
        try:
            named = self._win.identity(handle)
        finally:
            self._win.close(handle)
        if named != row.identity:
            raise HistoryLockError('evidence path now names a different file: ' + str(row.path))
        self._counts['identity_checks'] += 1

    def _read_digest(self, row, return_bytes=False):
        self._check_identity(row)
        if not self._win.seek(row.handle, 0, None, os.SEEK_SET):
            raise ctypes.WinError(ctypes.get_last_error())
        digest = hashlib.sha256()
        chunks = [] if return_bytes else None
        # Native handles avoid the CRT descriptor limit (commonly 512 and at
        # most 8192), so a session can hold more than 10,000 evidence files.
        buffer = self._read_buffer
        size = wintypes.DWORD()
        while True:
            if not self._win.read(row.handle, buffer, len(buffer), ctypes.byref(size), None):
                raise ctypes.WinError(ctypes.get_last_error())
            if not size.value:
                break
            chunk = ctypes.string_at(buffer, size.value)
            digest.update(chunk)
            self._counts['hash_bytes'] += len(chunk)
            if return_bytes:
                chunks.append(chunk)
        self._counts['hash_reads'] += 1
        self._check_identity(row)
        return digest.hexdigest(), b''.join(chunks) if return_bytes else None

    def _sha(self, value):
        self._ensure_open()
        self._counts['requests'] += 1
        path = self._path(value)
        key = self._key(path)
        if key in self._files:
            row = self._files[key]
            self._check_identity(row)
            self._counts['cache_hits'] += 1
            return row.digest
        self._lock_directories(path)
        handle = self._win.open(path)
        try:
            identity = self._win.identity(handle)
            row = _HeldFile(path, handle, identity, '')
            row.digest, _ = self._read_digest(row)
        except BaseException:
            self._win.close(handle)
            raise
        self._files[key] = row
        self._counts['unique_files'] += 1
        return row.digest

    def _abort(self):
        self._failed = True
        self.close()

    def sha(self, value):
        with self._mutex:
            try:
                return self._sha(value)
            except BaseException:
                self._abort()
                raise

    def verify(self, mapping):
        with self._mutex:
            try:
                self._ensure_open()
                if not hasattr(mapping, 'items'):
                    raise HistoryLockError('path to SHA256 mapping required')
                for value, digest in mapping.items():
                    path = self._path(value)
                    key = self._key(path)
                    if not isinstance(digest, str) or re.fullmatch(r'[0-9a-fA-F]{64}', digest) is None:
                        raise HistoryLockError('valid SHA256 expectation required: ' + str(value))
                    digest = digest.lower()
                    if key in self._expected and self._expected[key] != digest:
                        raise HistoryLockError('conflicting expected hashes: ' + str(path))
                    self._expected[key] = digest
                    if self._sha(path) != digest:
                        raise HistoryLockError('source/evidence hash mismatch: ' + str(path))
                return self.stats()
            except BaseException:
                self._abort()
                raise

    seal = verify

    def read_bytes(self, value):
        with self._mutex:
            try:
                digest = self._sha(value)
                row = self._files[self._key(self._path(value))]
                actual, payload = self._read_digest(row, return_bytes=True)
                if actual != digest:
                    raise HistoryLockError('evidence changed during read: ' + str(row.path))
                return payload
            except BaseException:
                self._abort()
                raise

    def final_audit(self):
        with self._mutex:
            try:
                self._ensure_open()
                for row in self._files.values():
                    digest, _ = self._read_digest(row)
                    if digest != row.digest:
                        raise HistoryLockError('final byte audit mismatch: ' + str(row.path))
                self._counts['final_audits'] += 1
                return self.stats()
            except BaseException:
                self._abort()
                raise

    def stats(self):
        with self._mutex:
            return dict(self._counts, locked_files=len(self._files),
                        locked_directories=len(self._directories),
                        expected_files=len(self._expected), closed=self._closed,
                        failed=self._failed, authority='LIVE_WINDOWS_DENY_WRITE_DELETE_HANDLES_V1')

    def close(self):
        with self._mutex:
            if self._closed:
                return
            self._closed = True
            for row in self._files.values():
                self._win.close(row.handle)
            self._files.clear()
            for handle, _ in reversed(list(self._directories.values())):
                self._win.close(handle)
            self._directories.clear()

    def __enter__(self):
        self._ensure_open()
        return self

    def __exit__(self, *_):
        self.close()

    def __del__(self):
        # CPython teardown may already have removed imported modules.
        try:
            self.close()
        except Exception:
            pass
