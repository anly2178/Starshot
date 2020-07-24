== in progress ==

TODO:
  LOTS OF TESTING - MATERIAL IS A NEW CLASS AND HASN'T BEEN TESTED MUCH
                  - this being said, it is based off of old working code; integration testing crucial here
  Add exceptions for invalid input w/ custom messages for the user
                  - e.g. invalid file loading (non-csv should raise exception)
  Add saving/loading functionality for Materials (saved as .txt is easiest)
                  - this NEEDS to be tested well; this is a large part of the usefulness of making a new Material class
                        - make following materials:
                              - gap
                              - SiO2
                              - GeO2
                              - Al2O3
                              - Si3N4
                          by loading data from .csv files
                  - EXAMPLE FILE SAVE:
                        SiO2_data.csv
                        gap_data.csv
                        
                        format: [material_name]_data.csv
                        
  Commenting + code cleanup
  Potentially add functionality for reading different file types?
