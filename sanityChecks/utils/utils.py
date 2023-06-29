"""
This is a collection of functions that don't really have a home
"""
import os


def get_project_root() -> str:
    current_directory = os.getcwd()
    files = os.listdir(current_directory)
    while "sanityChecks" not in files:
        current_directory = os.path.dirname(current_directory)
        files = os.listdir(current_directory)
    return current_directory
