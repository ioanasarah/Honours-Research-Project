version = "1.0.4"

''' Change logs
1.0.5:
        - Created new version of the GPR rules code which disambiguates all the possible OR rules and keeps all AND rules
                so that all AND rules can be evaluated at the same time and therefore can be trimmed
        - Added M-value based code for finding eligible genes that could be trimmed if trimming is selected for AND rules
        - New code does not deal with promiscuity and thus promiscuity leads to NotImplementedError 
        - Removed per_reaction_ datasets creation and saving
        - Removed creation of contribution datasets that resulted from promiscuity correction
        - Refactored transformation and thresholding code to its own function for better readabilityc
        - Added metafile saver for the specific settings used to create combinations folder
        
1.0.4:
        - Fixed bug where added metafile did not have .json extension and so would lead to an error when attempting .split(".")[1]
                in save_with_tries_function inside route_optimization.py
        -- Author: Jelle Bonthuis
        -- Date: 2024-12-19

1.0.3:
        - Fixed problem with GPR rules and value tokenization not working correctly (sum.15 would be the OR rule instead of just SUM). 
        - Added saving of partial datasets as .feather file for reading and histogram plotting later
        - Fixed issue with trimmin overtrimming when amount of genes were too small (only in not-for-promiscuity, 
                promiscuity related logic needs to be still fixed) (this only happened in the evaluate contribution but not evaluate
                expression
        -- Author: Jelle Bonthuis
        -- Date: 2024-12-10

1.0.2:  
        - Added a fix for a bug in the creation of taskstructure due to some numerals in excel being of int64 type.
        - Ensured that unsuccesful json file creation deletes the file.
        
        -- Author: Jelle Bonthuis
        -- Date: 2024-12-04
'''