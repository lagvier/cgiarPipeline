#------------------------------------------------------------------------------
# Name for the module (to identify in the modifications table): dcRaw
# 
# Functions for data consistency
# 1. Check names
# 2. Check consistency
#    2.1. Only informative (list of warnings): comprehensive
#    2.2. Modifications: only for yields but potentially extended to other traits
#         To be implemented later
# 3. Compute variables? (maybe in some other module, after loading data)
# 
# Output options for check consistency:
# - Only informative (option 2.1):
#   - list of warnings (ready)
#   - and visualization (would need modifications or new functions, to be implemented later)
# - Data modifications (option 2.2) => Tag inconsistencies and then change values to NA or 0
#                                      To be implemented later
# 
# Scope: Only inconsistencies as outliers are already detected in the QA module
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# 1. Check names
# 1.1. Identify names
#      - Convert CO numbers to short labels
#      - Output: text listing all traits not recognized
# 1.2. Renaming (matching old with new names).
#      Proposal: Use the menu Other functions / Trait Transformations
# 
# Process 1.2 will modify the data.
# Options:
# - Add new traits (columns) with the new names
# - Add matching list in the metadata table
#------------------------------------------------------------------------------

checkNames <- function(phenoDTfile, # The data object structure produced from bioflow
                       crop) {
  
  # Get internal copy of data and convert to short labels
  
  mydata <- phenoDTfile$data$pheno
  
  if (crop == 'potato')
    mydata <- suppressWarnings(st4gi::convert.co.pt(mydata, 'co.to.labels'))
  
  if (crop == 'sweetpotato')
    mydata <- suppressWarnings(st4gi::convert.co.sp(mydata, 'co.to.labels'))
  
  # Check names
  
  if (crop == 'potato')
    tmp <- st4gi::check.names.pt(mydata)
  
  if (crop == 'sweetpotato')
    tmp <- st4gi::check.names.sp(mydata)
  
  # Output: check.names.xx will produce a warning with a list of names not recognized
  # so the user can match those with standard labels
  
  return(tmp)

}

matchNames <- function(phenoDTfile, crop) {
  
  # This should give the option to match column names for traits
  # This should be reactive? produce changes in the checkNames function output?
  
}

#------------------------------------------------------------------------------
# 2. Check consistency
#    - This will spot potential problems
#    - Only inconsistencies will be reported as outliers' detection is in the QA module
#    Outputs in the Dashboard:
#    - List of inconsistencies in data.frame format.
#    - Matrix representation of the data.frame with colors to spot potential problems:
#      - Red cells: inconsistencies (impossible values)
#      This will need modification of the functions to point out specific cells
#------------------------------------------------------------------------------

checkData <- function(phenoDTfile, # The data.frame output from checkNames
                      crop, format = 'data.frame') {
  
  # Get internal copy of data after running checkNames and matchNames
  
  mydata <- phenoDTfile
  
  # Run check consistency
  # Use high number of f to avoid detecting outliers
  
  if (crop == 'potato')
    output <- suppressWarnings(st4gi::check.data.pt(mydata, f = 100, format))
  
  if (crop == 'sweetpotato')
    output <- suppressWarnings(st4gi::check.data.sp(mydata, f = 100, format))
  
  # Output could be text or data.frame

  # print(output) # Text output
  return(output) # Data frame output
  
}
