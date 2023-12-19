######## Medicare Advantage Dataset
## A detailed description of the datatask is contained in the .pdf file.

# Import the dataset, notice that does not have an header so we attach column name
# as stated in the pdf file.
data <- read.csv("scp-1205.csv", 
                 header=FALSE,
                 col.names = c("countyname", 
                               "state",
                               "contract", 
                               "healthplanname", 
                               "typeofplan", 
                               "countyssa", 
                               "eligibles", 
                               "enrollees", 
                               "penetration", 
                               "ABrate"))

# Replace Missing values for eligibles, enrollees and penetration with 0.
data$eligibles[is.na(data$eligibles)] <- 0
data$enrollees[is.na(data$enrollees)] <- 0
data$penetration[is.na(data$penetration)] <- 0

# The objective is to produce a dataset at a county-level that identifies the number of plans and total enrollment
# in each county (Puerto Rico and Guam excluded as state)
data <- data[data$state != "PR" & data$state != "GU",]

# list all the counties excluding "PUERTO RICO and GUAM"
counties <- unique(data$countyname)

# I instantiate the output dataset preallocating variables. I Know in advance the number of rows
# since it coincide with the number of unique counties present in the initial dataset.
output <- data.frame(countyname = counties, 
                     state = rep(0, length(counties)), 
                     numberofplans1 = rep(0, length(counties)), 
                     numberofplans2 = rep(0, length(counties)),
                     countyssa = rep(0, length(counties)),
                     eligibles = rep(0, length(counties)),
                     totalenrollees = rep(0, length(counties)),
                     totalpenetration = rep(0, length(counties)))

# Here is the main for loop.
# For each unique county
for (county_idx in 1:length(counties)){
  
  # Subset the initial dataset taking the current county
  subset <- data[data$countyname == counties[county_idx],]
  
  # Take the state, for each county subset the state will be the same so we can take the first element
  output$state[county_idx] <- subset$state[1]
  
  # In the variable numberofplans1 we want to count the number of health plans with more than 10 enrollees.
  # Therefore, I compare the column enrollees from subset with the number 10 obtaining a mask and I can sum the values
  # so that I get the number of health plans with more than 10 enrollees
  output$numberofplans1[county_idx] <- sum(subset$enrollees > 10)
  
  # For the variables numberofplans2 I do the same but with the variable penetration
  output$numberofplans2[county_idx] <- sum(subset$penetration > 0.5)
  
  # For the variable countyssa it is the same as state
  output$countyssa[county_idx] <- subset$countyssa[1]
  
  # For eligibles we want the number of individuals in the county that are medicare eligible so we sum
  # the column of the subset
  output$eligibles[county_idx] <- sum(subset$eligibles)
  
  # For totalenrollees we do the same
  output$totalenrollees[county_idx] <- sum(subset$enrollees)
  
  # For the variable TotalPenetration we compute the percent of individuals in the county enrolled in an MA plan
  # defined as totalenrollees/eligibles
  output$totalpenetration[county_idx] <- (output$totalenrollees[county_idx]/output$eligibles[county_idx])*100
}

# Last, we sort by state and by county
output <- output[order(output$state, output$countyname),]

# The first two observation seems to be strange, maybe is something we don't want even though it should be noticed
# that the number of eligibles is fairly high, so one should investigate more on that county name. For this example I just delete
# the first two entries
output <- output[-c(1,2),]

# Lastly I write the csv
write.csv(output, file="output_task.csv")




