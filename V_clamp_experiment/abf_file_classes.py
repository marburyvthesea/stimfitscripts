import sys
sys.path.append('/Users/johnmarshall/Documents/Analysis/PythonAnalysisScripts/jupyterscripts/');
import compile_by_abf_file

class analyzed_abf_file_pp40(object):
    ##class to store information for a stimfit .csv file
    
    def __init__(self, abf_file_name):
    #initialize file as name of able file
        self.abf_file_name = abf_file_name ; 
        
    def directory(self, root_dir):
    #returns the directory or directories containing the analyzed .csv file
        found_dir = compile_by_abf_file.find_dir_iterative_for_extension(root_dir, self.abf_file_name)[0]; 
        return(found_dir)

class V_clamp_experiment(object):
    ## class to store data on voltage clamp/KA application experiments
    #object is a list of files with 1st file being the baseline
    
    def __init__(self, expt_name, expt_type, raw_files):
        #initialize object as a list of the raw data files
        self.raw_files = raw_files; 
        self.expt_name = expt_name;
        self.expt_type = expt_type; 
        
    def analyze_raw(self, raw_files, search_dir):
        ##here call     
        #cycle through raw_files, create analyzed file class to associate with directory
        files_to_analyze = [];
        for fname in raw_files:
            stimfit_file = analyzed_abf_file_pp40(fname);
            files_to_analyze.append(str(stimfit_file.directory(search_dir))+'/'+str(fname)+'_stimfitanalysis.csv');
        print(files_to_analyze);
        analyze_experiment_from_csv.analyze_experiment_from_csv(files_to_analyze[0], str(self.expt_name), files_to_analyze[1:])
    
    def add_to_master(self, master_path):
        ##updates master df with information by loading master, updating, and resaving
        #does not save
        path_to_master_df = '/Users/johnmarshall/Documents/Analysis/RecordingData/masterabfspreadsheet.xlsx';
        master_df_ = pd.read_excel(path_to_master);
        for abf_file in self.raw_files: 
            stimfit_file = analyzed_abf_file_pp40(abf_file); 
            d = {'associated experiment' : pd.Series([self.expt_name], index=['a']),
            'directory' : pd.Series([str(stimfit_file.directory(root_dir))], index=['a']),
            'experiment type': pd.Series([self.expt_type], index=['a']),
            'file': pd.Series([stimfit_file.abf_file_name], index=['a'])
            };
            to_add = pd.DataFrame(d);
            master_df_ = pd.concat([master_df_, to_add], axis = 0);
        return(master_df_)
    
    def update_master_file(self, path_to_master, path_to_save):
        ##updates and saves master df 
        updated_master_df = self.add_to_master(path_to_master); 
        updated_master_df.to_excel(path_to_save); 
        return(); 