
function mustBeFileType(file_path,needed_file_format_ending)
    arguments
        file_path (1,1) string
        needed_file_format_ending (1,1) string
    end
            assert(exist(file_path,'file') ==2 , "Does the file exist ? Check again!")
            assert(~isempty(regexp(file_path,needed_file_format_ending + "$")),...
                   "Input must be a " + needed_file_format_ending + " file!!")
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
