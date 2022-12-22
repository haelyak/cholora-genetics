import unittest
import cholera as hw

class SimpleTests(unittest.TestCase):

    def testA(self):
        # hasHomolog
        self.assertTrue(hw.hasHomolog(m1,m2,3)==[1, 1, 1],msg="hasHomolog(m1,m2,3) failed.")

    def testB(self):
        # hasHomolog
        self.assertTrue(hw.hasHomolog(m1,m2,5)==[1, 0, 1],msg="hasHomolog(m1,m2,5) failed.")

    def testC(self):
        # hasHomolog
        self.assertTrue(hw.hasHomolog(m1,m2,8)==[0, 0, 1],msg="hasHomolog(m1,m2,8) failed.")

    def testD(self):
        # hasHomolog
        self.assertTrue(hw.hasHomolog(m1,m2,10)==[0, 0, 1],msg="hasHomolog(m1,m2,10) failed.")


    def testE(self):
        # hasHomolog
        self.assertTrue(hw.hasHomolog(m1,m2,11)==[0, 0, 0],msg="hasHomolog(m1,m2,11) failed.")

    def testF(self):
        # islands
        self.assertTrue(hw.islands([1,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1],4)==[(14, 21), (8, 12)],msg="islands([1,0,0,0,1,0,1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1],4) failed.")


    def testG(self):
        # islands
        self.assertTrue(hw.islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],5)==[(13, 20), (0, 5)],msg="islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],5) failed.")

        
    def testH(self):
        # islands
        self.assertTrue(hw.islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],6)==[(13, 20)],msg="islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],6) failed.")


    def testI(self):
        # islands
        self.assertTrue(hw.islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],10)==[],msg="islands([0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0],10) failed.")

        
if __name__ == '__main__':

    m1=[[3,6,6,7],
        [2,2,4,3],
        [4,6,5,10]]

    m2=[[1,3,4,1],
        [2,4,3,1],
        [8,5,6,1]]

    unittest.main()
