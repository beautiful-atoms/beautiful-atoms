import pytest
import sys
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def get_pytest_args():
    """Extract pytest arguments from Blender's command-line arguments."""
    try:
        # Find where "--" appears in the args (separates Blender args from pytest args)
        index = sys.argv.index("--")
        return sys.argv[index + 1 :]  # Extract pytest arguments
    except ValueError:
        logging.warning(
            "No `--` separator found. Running pytest with default arguments."
        )
        return []


def main():
    """Run pytest inside Blender and exit with the correct status code."""
    pytest_args = get_pytest_args()

    logging.info(f"Running tests with arguments: {pytest_args}")

    # Run pytest with extracted arguments
    result = pytest.main(pytest_args)

    if result.value != 0:
        logging.error(f"Tests failed with exit code: {result.value}")
        sys.exit(result.value)
    else:
        logging.info("All tests passed successfully.")


if __name__ == "__main__":
    main()
