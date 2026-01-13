classdef medium
    
    
    properties
        file_names
        file_location
        medium_composition
        manual_set_boundaries
    end
    
    methods
        function obj = medium(filenames,filelocation)
            %arguments
            %       filenames (1,:) string {mustBeFileType(filenames,"tsv")}
            %       filelocation (1,:) string
            %end
            obj.file_names = filenames;
            
            if length(filenames) == length(filelocation)
                obj.file_location = filelocation;
            elseif length(filelocation) ==1
                obj.file_location = filelocation(ones(1,length(filenames)));
            else
                error("The filelocation arguments needs to be one string, or as many strings as files are given!")
            end     
        end
        
        function obj = read_medium_files(obj, column_to_join_by,columns_in_both)
           arguments
                 obj
                 column_to_join_by
                 columns_in_both =["Mets_Recon3D","ExRxns_Recon3D","ExRxns_HumanGEM","Mets_HumanGEM"]
           end
           col_postfix = regexprep(obj.file_names,".tsv","");
           column_to_join = column_to_join_by + "_" + col_postfix;
           med = {};
           for idx_file=1:length(obj.file_names)
               file = obj.file_location(idx_file) + filesep + obj.file_names(idx_file);
               med{idx_file} = readtable(file, ...
                                         'FileType', 'text', 'Delimiter', '\t');               
               med{idx_file}.source = file(ones(size(med{idx_file},1),1))
               
               med{idx_file}.Properties.VariableNames(string(med{idx_file}.Properties.VariableNames) == column_to_join_by(idx_file)) = column_to_join(idx_file);
           end
           
           tab = med{1};
           col_conc = column_to_join(1);
           for med_idx=2:length(obj.file_names)
                medium = outerjoin(tab,med{med_idx},...
                                   'MergeKeys',true,'Keys',...
                                   columns_in_both);
                                    
                medium{find(isnan(medium{:, [col_conc]})),col_conc} = 0;
                medium{find(isnan(medium{:, [column_to_join(med_idx)]})),column_to_join(med_idx)} = 0;
                
                % adding up the metabolites, in case they are in both media components
                medium{:, ["Concentration_M"]} = medium{:, [col_conc]} + medium{:, [column_to_join(med_idx)]};
                medium{:, ["Concentration_mM"]} = medium{:, ["Concentration_M"]}*1000;
                medium.(col_conc) = [];
                medium.(column_to_join(med_idx)) = []; 
                col_conc = "Concentration_M_tab";
                medium.(col_conc) = medium.Concentration_M;
                tab = medium;
           end
           medium.(col_conc) = [];
           obj.medium_composition = medium;
            
        end
        
        function obj = add_additional_rxns_boundaries(obj, unwanted_import, unwanted_export,wanted_import, wanted_export)
           arguments
              obj
              unwanted_import (1,:) string =[]
              unwanted_export (1,:) string =[]
              wanted_import (1,:) string =[]
              wanted_export (1,:) string =[]
           end
             obj.manual_set_boundaries.unwanted_import = unwanted_import;
             obj.manual_set_boundaries.unwanted_export = unwanted_export;
             obj.manual_set_boundaries.wanted_import = wanted_import;
             obj.manual_set_boundaries.wanted_export = wanted_export;
           
        end
    end
end

function mustBeFileType(file_path,needed_file_format_ending)
    arguments
        file_path (1,1) string
        needed_file_format_ending (1,1) string
    end
            assert(exist(file_path,'file') ==2 , "Does the file exist ? Check again!")
            assert(~isempty(regexp(file_path,needed_file_format_ending + "$")),...
                   "Input must be a " + needed_file_format_ending + " file!!")
end


