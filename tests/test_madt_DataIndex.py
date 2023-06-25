import unittest
from madt.madt_DataIndex import mkdir 

class TestDataIndex(unittest.TestCase):

	def test_data_index(self):
		self.assertEqual(mkdir('test'), 'test')

unittest.main()