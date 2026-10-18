"""Old suites and retained FAIL plus bounded-refinement strategy checks."""
from run_process_event_refined import configure, EXTRA_SOURCES
import check_process_event as checker


def main():
    configure()
    checker.TESTS = checker.TESTS + ['event_response_refined', 'process_event_refined_runner']
    checker.REQUIRED = checker.REQUIRED | set(EXTRA_SOURCES)
    return checker.main()


if __name__ == '__main__':
    raise SystemExit(main())
