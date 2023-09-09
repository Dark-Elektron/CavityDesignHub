import unittest
from utils.shared_functions import *


class TestSharedFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls) -> None:
        # load shape space
        with open(r"D:\Dropbox\CavityDesignHub\unittest_files\reference_files\test_shape_space.json", 'r') as file:
            cls.shape_space = json.load(file)['C3795']
        with open(r"D:\Dropbox\CavityDesignHub\unittest_files\reference_files\test_shape_space_scale_2.json") as file:
            cls.shape_space_scale_2 = json.load(file)['C3795']

    @classmethod
    def tearDownClass(cls) -> None:
        pass

    def setUp(self) -> None:
        pass

    def tearDown(self) -> None:
        pass

    def test_f2b_slashes(self):
        if os.name == 'nt':
            self.assertEqual(f2b_slashes('D:/User/directory'), 'D:\\User\\directory')
        else:
            self.assertEqual(f2b_slashes('D:\\User\\directory'), 'D:/User/directory')

    def test_scale_cavity_geometry(self):
        self.assertEqual(scale_cavity_geometry(self.shape_space, 2), self.shape_space_scale_2)


if __name__ == '__main__':
    unittest.main()
