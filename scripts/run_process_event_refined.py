"""Named bounded force-difference refinement; original physical driver intact.

Install a numerical response strategy only. The old driver remains reproducible
and its failed cloud is never continued or relabelled. No physical map is patched.
"""
import run_process_event as driver
from afsd_backend.pysph.event_response_refined import response_from_trace

EXTRA_SOURCES = [
    'scripts/run_process_event_refined.py', 'scripts/check_process_event_refined.py',
    'src/afsd_backend/pysph/event_response_refined.py',
    'tests/review/test_event_response_refined.py',
    'tests/review/test_process_event_refined_runner.py',
    'docs/PROCESS_EVENT_REFINEMENT_CONTRACT.md',
    'docs/PROCESS_EVENT_REFINED_ROOT_CONTRACT.md',
]


def configure():
    driver.response_from_trace = response_from_trace
    driver.SOURCES = list(dict.fromkeys(driver.SOURCES + EXTRA_SOURCES))
    driver.SETTINGS = dict(driver.SETTINGS,
        preconditioner='FROZEN_EVENT_B_C_SCHUR_BOUNDED_FORCE_DIFFERENCE_REFINEMENT',
        force_difference_maximum_halvings=2)
    return driver


if __name__ == '__main__':
    raise SystemExit(configure().main())
