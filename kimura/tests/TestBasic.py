import unittest, os, shutil 
import kimura

#Original input file 
TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'oocyte.txt')
OUTPUT_TEMP_DIR = os.path.join(os.path.dirname(__file__), 'temp/')

class TestBasic(unittest.TestCase):
    try:
        os.mkdir(OUTPUT_TEMP_DIR)
    except OSError as err:
        print("OS error: {0}".format(err))
    else:
        print ("Successfully created the directory %s " % OUTPUT_TEMP_DIR)
    p = kimura.run(TESTDATA_FILENAME, OUTPUT_TEMP_DIR + 'test', 10, 10, 100)

    #set number of simulated dataset in each run 1000, p-value can't be below 0.001
    def test_p_value(self):
        self.assertTrue(self.p >= 0.001)

    def test_check_output_files(self):
        self.assertTrue(os.path.isfile(OUTPUT_TEMP_DIR + 'test.kimura.comp.csv'))
        self.assertTrue(os.path.isfile(OUTPUT_TEMP_DIR + 'test.kimura.simu.csv'))
        self.assertTrue(os.path.isfile(OUTPUT_TEMP_DIR + 'test.monte_carlo.tsv'))
        self.assertTrue(os.path.isfile(OUTPUT_TEMP_DIR + 'test.observed.comp.csv'))
        self.assertTrue(os.path.isfile(OUTPUT_TEMP_DIR + 'test.summary.txt'))

if __name__ == '__main__':
    unittest.main()
