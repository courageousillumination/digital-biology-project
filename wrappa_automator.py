#!/usr/bin/env python
"""
This is a script that will run wrappa (through the web interface) on a large
set of files.
"""

import logging
import os

from pdb_parser import find_pdb_files
from selenium import webdriver


def save_page_as(browser, file_name):
    """
    Save the output of a page.
    """

    with open(file_name, "w") as fout:
        fout.write(browser.find_element_by_tag_name("pre").text)

def run_wrappa(browser, pdb_file):
    """
    This will run through all of the web interface using selenium. The input
    should be a pdb_file full path. It will download the wrappers/bonds as
    PDB_NAME_wrappers.txt and PDB_NAME_bonds.txt.
    """

    # Wrappa has 3 MB limit
    if os.path.getsize(pdb_file) > 3000000:
        logging.warn("%s is too large (size is %d), skipping", pdb_file, os.path.getsize(pdb_file))
        return False
    if os.path.isfile(pdb_file[:-4] + "_bonds.txt"):
        logging.warn("%s has already been processed, skipping", pdb_file)
        return False

    full_path = os.path.abspath(pdb_file)
    directory, file_name = os.path.split(full_path)
    pdb_name = file_name[:-4]

    try:
        browser.get("http://www.wrappa.org/wrappa01/wrappa")
        browser.find_element_by_name("pdbFileName").send_keys(full_path)
        browser.find_element_by_xpath("//*[@type='submit']").click()
        # Use the default configuration
        browser.find_element_by_xpath("//*[@type='submit']").click()
        # Analyze
        browser.find_element_by_xpath("//*[@type='submit']").click()
        # Download the files created
        browser.find_element_by_link_text("Bonds").click()
        save_page_as(browser, os.path.join(directory, pdb_name + "_wrappers.txt"))
        browser.back()
        browser.find_element_by_link_text("Wrappers").click()
        save_page_as(browser, os.path.join(directory, pdb_name + "_bonds.txt"))
    except Exception as e:
        logging.exception(e)
        logging.error("Encounterd exception for %s", pdb_file)
        return False

    return True

def process_directory(browser, directory):
    """
    Process an entire directory.
    """

    count = 0
    for pdb in find_pdb_files(directory):
        if run_wrappa(browser, pdb):
            logging.info("Processed %s", pdb)
            count += 1
    logging.info("Fully processed %d pdb files", count)

def initialize():
    """
    Initialize the browser and logging
    """

    logging.basicConfig(level=logging.INFO)

    browser = webdriver.Firefox()
    browser.get("http://www.wrappa.org/")
    browser.find_element_by_name("termsAccepted").click()
    browser.find_element_by_name("termsAccepted").submit()
    return browser

def main():
    """
    Run the main functionality of this script.
    """

    browser = initialize()
    process_directory(browser, "data")
    browser.close()

if __name__ == "__main__":
    main()
