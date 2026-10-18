"""Explicit FP64 CUDA linear algebra; no CPU fallback or global monkeypatch.

PySPH, particle storage, matrix assembly and posterior physics gates stay on
CPU. Calls include both transfers and synchronize before reporting wall time.
Dependencies are isolated from the frozen compiled PySPH environment.
"""
from pathlib import Path
import importlib.metadata
import json
import os
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OVERLAY = ROOT/'.envs/gpu_overlay_20261007'
_cp = None
_dll_handles = []
_calls = []


def initialize():
    global _cp
    if _cp is None:
        import ctypes
        if not OVERLAY.is_dir():
            raise RuntimeError('Project CUDA dependency overlay is missing')
        sys.path.insert(0, str(OVERLAY))
        # The desktop inherits a different Anaconda base containing CUDA12.1.
        # Bind discovery to THIS interpreter, not that unrelated environment.
        os.environ['CONDA_PREFIX'] = sys.prefix
        os.environ.setdefault('CUPY_CACHE_DIR', str(ROOT/'outputs/gpu_acceleration_20261007/kernel_cache'))
        toolkit = Path(os.environ.get('AFSD_CUDA_TOOLKIT', 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.8')).resolve()
        if os.name != 'nt' or not (toolkit/'bin').is_dir():
            raise RuntimeError('This explicit CUDA DLL binding requires the registered Windows toolkit')
        # PATH/add_dll_directory alone are insufficient: cuda-pathfinder can
        # discover the inherited Anaconda 12.1 DLLs before CUDA_PATH. Bind the
        # actual libraries before importing CuPy and refuse a polluted process.
        os.environ['CUDA_PATH'] = os.environ['CUDA_HOME'] = str(toolkit)
        os.environ['PATH'] = str(toolkit/'bin')+os.pathsep+os.environ.get('PATH', '')
        _dll_handles.append(os.add_dll_directory(str(toolkit/'bin')))
        kernel32 = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel32.GetModuleHandleW.argtypes = [ctypes.c_wchar_p]
        kernel32.GetModuleHandleW.restype = ctypes.c_void_p
        kernel32.GetModuleFileNameW.argtypes = [ctypes.c_void_p, ctypes.c_wchar_p, ctypes.c_uint32]
        kernel32.GetModuleFileNameW.restype = ctypes.c_uint32

        def module_path(handle):
            buffer = ctypes.create_unicode_buffer(32768)
            length = kernel32.GetModuleFileNameW(handle, buffer, len(buffer))
            if not length or length >= len(buffer):
                raise RuntimeError('Cannot verify loaded CUDA DLL path')
            return Path(buffer.value).resolve()

        loaded_paths = {}
        for name in ('cudart64_12.dll', 'nvJitLink_120_0.dll', 'cublasLt64_12.dll',
                     'cublas64_12.dll', 'cusparse64_12.dll', 'cusolver64_11.dll',
                     'nvrtc-builtins64_128.dll', 'nvrtc64_120_0.dll'):
            expected = toolkit/'bin'/name
            if not expected.is_file():
                raise RuntimeError('Registered CUDA DLL missing: '+str(expected))
            existing = kernel32.GetModuleHandleW(name)
            if existing and module_path(existing) != expected:
                raise RuntimeError('Conflicting CUDA DLL already loaded; use a fresh process: '+str(module_path(existing)))
            # DLL_LOAD_DIR | SYSTEM32: use the pinned DLL's dependencies and
            # Windows system libraries, not unrelated Anaconda DLL directories.
            handle = ctypes.WinDLL(str(expected), winmode=0x00000100 | 0x00000800)
            _dll_handles.append(handle)
            actual = module_path(handle._handle)
            if actual != expected:
                raise RuntimeError('CUDA DLL origin mismatch: '+str(actual))
            loaded_paths[name] = str(actual)
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError('CUDA GPU unavailable; CPU fallback is forbidden')
        cp.cuda.Device(0).use()
        if tuple(cp.cuda.nvrtc.getVersion()) < (12, 8):
            raise RuntimeError('This Blackwell GPU requires the registered NVRTC 12.8 or newer')
        initialize._binding = dict(toolkit=str(toolkit), loaded_dll_paths=loaded_paths,
            process_environment_only=True, conflicting_loaded_libraries_rejected=True)
        _cp = cp
    return _cp


def _input(matrix, rhs=None):
    if not isinstance(matrix, np.ndarray) or matrix.dtype != np.float64:
        raise ValueError('FP64 NumPy matrix required; no precision conversion')
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not matrix.size or not np.isfinite(matrix).all():
        raise ValueError('finite nonempty square matrix required')
    if rhs is not None:
        if not isinstance(rhs, np.ndarray) or rhs.dtype != np.float64 or rhs.ndim not in (1, 2):
            raise ValueError('FP64 NumPy vector/matrix RHS required')
        if rhs.shape[0] != matrix.shape[0] or not np.isfinite(rhs).all():
            raise ValueError('finite conforming RHS required')


def _run(kind, matrix, rhs=None):
    _input(matrix, rhs)
    cp = initialize()
    began = time.perf_counter()
    try:
        import cupyx
        with cupyx.errstate(linalg='raise'):
            device_matrix = cp.asarray(matrix)
            answer = (cp.linalg.inv(device_matrix) if kind == 'inv' else
                      cp.linalg.solve(device_matrix, cp.asarray(rhs)))
            result = cp.asnumpy(answer)
        cp.cuda.get_current_stream().synchronize()
        if result.dtype != np.float64 or not np.isfinite(result).all():
            raise np.linalg.LinAlgError('nonfinite or non-FP64 CUDA answer')
    except Exception as error:
        _calls.append(dict(operation=kind, n=matrix.shape[0], wall_s=time.perf_counter()-began, passed=False))
        raise np.linalg.LinAlgError('Explicit CUDA '+kind+' failed: '+repr(error)) from error
    _calls.append(dict(operation=kind, n=matrix.shape[0], wall_s=time.perf_counter()-began, passed=True))
    return result


def solve(matrix, rhs):
    return _run('solve', matrix, rhs)


def inv(matrix):
    return _run('inv', matrix)


def statistics(reset=False):
    rows = list(_calls)
    if reset:
        _calls.clear()
    groups = {}
    for row in rows:
        key = row['operation']+'_'+str(row['n'])
        group = groups.setdefault(key, dict(calls=0, wall_s=0., failures=0))
        group['calls'] += 1
        group['wall_s'] += row['wall_s']
        group['failures'] += int(not row['passed'])
    return dict(backend='CuPy CUDA FP64', cpu_fallback=False, calls=len(rows), groups=groups)


def environment():
    cp = initialize()
    from cuda.pathfinder import load_nvidia_dynamic_lib
    resolved = {}
    for name in ('cudart', 'nvrtc', 'cusolver', 'cublas', 'cusparse'):
        actual = Path(load_nvidia_dynamic_lib(name).abs_path).resolve()
        if str(actual) not in initialize._binding['loaded_dll_paths'].values():
            raise RuntimeError('CuPy resolved an unregistered CUDA DLL: '+str(actual))
        resolved[name] = str(actual)
    props = cp.cuda.runtime.getDeviceProperties(0)
    name = props['name']
    return dict(backend='CuPy CUDA FP64', device=0,
        device_name=name.decode() if isinstance(name, bytes) else str(name),
        total_memory_bytes=int(props['totalGlobalMem']),
        compute_capability=cp.cuda.Device(0).compute_capability,
        driver_version=cp.cuda.runtime.driverGetVersion(), runtime_version=cp.cuda.runtime.runtimeGetVersion(),
        cupy_version=cp.__version__, numpy_version=np.__version__, python=sys.version,
        nvrtc_version=cp.cuda.nvrtc.getVersion(), cuda_toolkit=os.environ.get('CUDA_PATH'),
        cuda_dll_binding=initialize._binding, cupy_resolved_dll_paths=resolved,
        overlay=str(OVERLAY), precision='float64', cpu_fallback=False)


def component_checks():
    rng = np.random.default_rng(20261007)
    checks = []
    for n in (3, 32, 97):
        a = rng.normal(size=(n, n)); a += (2*n)*np.eye(n)
        b = rng.normal(size=n)
        cpu = np.linalg.solve(a, b); gpu = solve(a, b)
        error = float(np.max(np.abs(cpu-gpu)))
        residual = float(np.linalg.norm(a@gpu-b, np.inf)/max(1., np.linalg.norm(b, np.inf)))
        inverse = inv(a)
        inverse_error = float(np.linalg.norm(a@inverse-np.eye(n), np.inf))
        assert error <= 1e-12 and residual <= 1e-12 and inverse_error <= 1e-12
        checks.append(dict(n=n, solution_abs_error=error, relative_residual=residual, inverse_residual=inverse_error))
    for a, b in ((np.zeros((3, 3)), np.ones(3)), (np.eye(3, dtype=np.float32), np.ones(3))):
        try:
            solve(a, b)
        except (ValueError, np.linalg.LinAlgError):
            pass
        else:
            raise AssertionError('singular or lower precision input must fail')
    return dict(status='PASS_CUDA_FP64_COMPONENTS', checks=checks,
                singular_and_float32_rejected=True, environment=environment(), statistics=statistics())


if __name__ == '__main__':
    report = component_checks()
    if len(sys.argv) != 2:
        raise SystemExit('explicit new report path required')
    target = Path(sys.argv[1])
    with target.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(report))
