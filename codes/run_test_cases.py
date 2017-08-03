#!/usr/bin/env python

import os   
from master_run import Master

class Run_Test(object):
   
    """THIS CODE USES THE MASTER CODE TO CREATE DEFAULT FILES WHICH ARE THEN
    COMPARED AGAINST PREVIOUSLY WRITTEN FILES (KNOWN TO WORK.)
    """
    
    def __init__(self, test_case='11fe'):
        
        self.test_case = test_case
        
        self._created_ymlfiles_list = None
        self._yml_orig, self._yml_new = None, None
        self._abun_orig, self._abun_new = None, None
        self._dens_orig, self._dens_new = None, None

        os.system('clear')
        print '\n\n\n'
        print '****************************************************'
        print '************ RUNNING TEST CASE FOR ' + self.test_case +' ************'         
        print '****************************************************'
        print '\n'
    
        self.run_test()
    
    def make_files(self):
        
        self._created_ymlfiles_list = Master(
          event=self.test_case, case='test', StoN='high',
          flag_run_simulation=False, run_uncertainties=False,
          flag_compute_features=False, make_kromer=False,
          flag_display_interface=False, verbose=False).run_master() 
        
    def get_dirs(self):
        
        if self.test_case == '11fe':
            top_dir_orig = os.path.abspath('./../test_cases/11fe') + '/'
            self._yml_orig = top_dir_orig + 'loglum-9.544.yml'
            self._abun_orig = (top_dir_orig + 'abundance_19.1_day.dat')
            self._dens_orig = (top_dir_orig + 'density_es-1.0_ms-1.0.dat')
                              
        elif self.test_case == '05bl':
            top_dir_orig = os.path.abspath('./../test_cases/05bl') + '/'
            self._yml_orig = (top_dir_orig + 'velocity_start-8100_loglum-8.617'
                              + '_time_explosion-12.0.yml')
            self._abun_orig = (top_dir_orig + 'abundance_es-0.7_ms-1.0_12.0'
                               + '_day.dat')
            self._dens_orig = (top_dir_orig + 'density_es-0.7_ms-1.0_12.0'
                               + '_day.dat')                              
                          
        self._yml_new = self._created_ymlfiles_list[0]
        path_new = os.path.dirname(self._yml_new)
        for fname in os.listdir(path_new):
            if fname[0:9] == 'abundance':
                self._abun_new = path_new + '/' + fname
            elif fname[0:7] == 'density':
                self._dens_new = path_new + '/' + fname        
    
    def compare_files(self):
        
        for fname, f1, f2 in zip(
                            ['ABUNDANCE', 'DENSITY', 'YML'],
                            [self._abun_orig, self._dens_orig, self._yml_orig],
                            [self._abun_new, self._dens_new, self._yml_new]
                            ):
        
            print '\n\n---->COMPARING ' + fname + ' FILES:\n\n'
            with open(f1, 'r') as inp1, open(f2, 'r') as inp2:
                for i, (line1, line2) in enumerate(zip(inp1, inp2)):
                    if line1 != line2:
                        print '-------->LINE ' + str(i + 1) +':\n'
                        print '    |---->OLD', line1.strip('\n') 
                        print '    |---->NEW', line2.strip('\n')
                        print '\n' 
    
        
    def run_test(self):
        self.make_files()
        self.get_dirs()
        self.compare_files()
                    

#master_obj = Run_Test(test_case='11fe')      
master_obj = Run_Test(test_case='05bl')      
