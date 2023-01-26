#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 12:29:24 2022

@author: prowe
"""


from antarc.pyr_save_stand_files import pyr_save_stand_files


class TestPyrSaveStandFiles:
    """
    Test that files are created by pyr_save_stand_files
    """

    def test_swd_v1_2017(self):
        """
        Test that files are created for happy SWD case
        """
        paramfile = "tests.parameters.swd_orig_params"
        pyr_save_stand_files(paramfile, "v1", 2017)

    def test_swd_v1_2015(self):
        """
        Test that files are not created for SWD when year is not in range
        """
        paramfile = "tests.parameters.swd_orig_params"

        # with pytest.raises(FileNotFoundError):
        #    assert pyr_save_stand_files(paramfile, "v1", 2015)
        pyr_save_stand_files(paramfile, "v1", 2015)

    def test_swd_v2_2022(self):
        """
        Test that files are not created for SWD when year is not in range
        """
        paramfile = "tests.parameters.swd_orig_params"

        # with pytest.raises(FileNotFoundError):
        #    assert pyr_save_stand_files(paramfile, "v1", 2015)
        pyr_save_stand_files(paramfile, "v2", 2022)

    def test_lwd_v1(self):
        """
        Test LWD files for year in range
        """
        paramfile = "tests.parameters.lwd_orig_params"
        pyr_save_stand_files(paramfile, "v1", 2020)

    def test_lwd_v2(self):
        paramfile = "tests.parameters.lwd_orig_params"
        pyr_save_stand_files(paramfile, "v2", 2022)


if __name__ == "__main__":
    testPyrSaveStandFiles = TestPyrSaveStandFiles()
    testPyrSaveStandFiles.test_lwd_v1()
