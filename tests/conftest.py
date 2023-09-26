import pytest
import shutil
import os


# Define a custom fixture to clean up directories before starting tests
@pytest.fixture(scope="session", autouse=True)
def clean_directories():
    # Define the directories to remove
    directories_to_remove = [
        "tests/pipeline-1",
        "tests/pipeline-2",
        # Add more directories to remove as needed
    ]

    # Remove the specified directories
    for directory in directories_to_remove:
        if os.path.exists(directory):
            shutil.rmtree(directory, ignore_errors=True)
