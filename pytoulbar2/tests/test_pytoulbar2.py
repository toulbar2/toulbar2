from unittest import TestCase
import pytoulbar2

class TestExtension(TestCase):
   def test_1(self):
      myCFN = pytoulbar2.CFN(2)
      sol = myCFN.Solve()
      self.assertEqual(sol[0],[])
      self.assertEqual(sol[1],0.0)
      self.assertEqual(sol[2],1)
 
