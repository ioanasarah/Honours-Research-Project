library(openxlsx)
library(stringr)

#PARAMETER 
metPrefix <- "MAM" #prefix of metabolites in the used model
outputFilePath <- "../../../../Data/MetabolicTasks_MACSBIO_v0.4_RAVEN.xlsx"#Path of the file to be written
#Also note that if the compartment is not indicated by a SINGLE letter after the ID, the scipt will fail

tasksCobra <- openxlsx::read.xlsx("../../../../Data/MetabolicTasks_MACSBIO_v0.4.xlsx",sheet = 1,sep.names = " ")
metIndex <- openxlsx::read.xlsx("../../../../Data/MetabolicTasks_MACSBIO_v0.4.xlsx",sheet = 2,sep.names = " ")
#cutting of unused cols
isNA <- which(colnames(tasksCobra) == "NA.")
if(length(isNA) != 0){
  tasksCobra <- tasksCobra[,-c(isNA:ncol(tasksCobra))]
}
#Checking the presence of the fields:
expectedCols <- c("ID","SYSTEM","SUBSYSTEM","DESCRIPTION","SHOULD FAIL","IN","IN LB","IN UB","OUT","OUT LB","OUT UB","EQU","EQU LB","EQU UB","COMP","COMMENTS")
if(sum(expectedCols == colnames(tasksCobra)) != ncol(tasksCobra)){
  print("THE STRUCTURE OF THE EXCEL SHEET IS NOT COMPATIBLE WITH THIS SCRIPT,")
  print("PLEASE MODIFY IT TO FOLLOW THE STRUCTURE LAYED OUT BY 'expectedCols' ")
}

#The task file should not be any bigger than the cobra version
tasksRaven <- as.data.frame(matrix(ncol = 17,nrow = nrow(tasksCobra)))
colnames(tasksRaven) <- c("ID","DESCRIPTION","SOULD FAIL","IN","IN LB","IN UB","OUT","OUT LB","OUT UB","EQU","EQU LB","EQU UB","CHANGED RXN","CHANGED LB","CHANGED UB","PRINT FLUX","COMMENTS")

#getting the position of each cobra task
taskID <- unique(tasksCobra$ID)
taskID <- taskID[!is.na(taskID)]
taskIndPos <- c(1:nrow(tasksCobra))[!is.na(tasksCobra$ID)]
taskIndEnd <- c(taskIndPos[-1]-1,nrow(tasksCobra))


#main loop
for (i in taskID){
  print(i)
  taskStart <- taskIndPos[i]
  taskEnd <- taskIndEnd[i]
  inputInds <- c(taskStart:taskEnd)[!is.na(tasksCobra$IN[taskStart:taskEnd])]
  #Input loop
  for (inInd in inputInds){
    inputID <- tasksCobra$IN[inInd]
    if(inputID == "ALLMETS"){
      tasksRaven[inInd,4] <- "ALLMETS"
    }else{
      comp <- substr(inputID,nchar(inputID)-1,nchar(inputID)-1)
      inputName <- metIndex[metIndex$Met_ID == inputID,1]
      
      RavenName <- paste0(inputName,"[",comp,"]")
      tasksRaven[inInd,4] <- RavenName #IN
    }
    
    tasksRaven[inInd,5] <- tasksCobra$`IN LB`[inInd] #IN LB
    tasksRaven[inInd,6] <- tasksCobra$`IN UB`[inInd] #IN UB
  }
  
  outputInds <- c(taskStart:taskEnd)[!is.na(tasksCobra$OUT[taskStart:taskEnd])]
  #Output loop
  for (outInd in outputInds){
    outputID <- tasksCobra$OUT[outInd]
    if(outputID == "ALLMETS"){
      tasksRaven[outInd,7] <- "ALLMETS" #OUT
    }else{
      comp <- substr(outputID,nchar(outputID)-1,nchar(outputID)-1)
      outputName <- metIndex[metIndex$Met_ID == outputID,1]
      
      RavenName <- paste0(outputName,"[",comp,"]")
      tasksRaven[outInd,7] <- RavenName #OUT
    }
    
    tasksRaven[outInd,8] <- tasksCobra$`OUT LB`[outInd] #OUT LB
    tasksRaven[outInd,9] <- tasksCobra$`OUT UB`[outInd] #OUT UB
  }
  
  #EQU
  hasEQU <- !is.na(tasksCobra$EQU[taskStart])
  if(hasEQU){
    
    equation <- tasksCobra$EQU[taskStart]
    leftHand <- gsub("=.*$","",equation)
    IDstarts <- str_locate_all(leftHand,metPrefix)
    IDstarts <- IDstarts[[1]]
    
    directionSign <- gsub("\\+","",paste0(str_extract_all(equation,"[:symbol:]")[[1]],collapse = ""))
    
    rightHand <- gsub("^.*=","",equation)
    IDstarts2 <- str_locate_all(rightHand,metPrefix)
    IDstarts2 <- IDstarts2[[1]]
    
    substrateNames <- apply(IDstarts,1, function(x){
      substrateID <- substr(leftHand,x[1],x[1]+8)
      comp <- substr(substrateID,nchar(substrateID)-1,nchar(substrateID)-1)
      paste0(metIndex$Met_Names[metIndex$Met_ID == substrateID],"[",comp,"]")
    })
    productNames <- apply(IDstarts2,1, function(x){
      productID <- substr(rightHand,x[1],x[1]+8)
      comp <- substr(productID,nchar(productID)-1,nchar(productID)-1)
      paste0(metIndex$Met_Names[metIndex$Met_ID == productID],"[",comp,"]")
    })
    
    equationNames <- paste0(paste0(substrateNames,collapse = " + ")," ",directionSign," ",paste0(productNames,collapse = " + "))
    tasksRaven[taskStart,10] <- equationNames #EQU
    tasksRaven[taskStart,11] <- tasksCobra$`EQU LB`[taskStart] #EQU LB
    tasksRaven[taskStart,12] <- tasksCobra$`EQU UB`[taskStart] #EQU UB
  }
  #General info
  tasksRaven[taskStart,1] <- i #ID
  tasksRaven[taskStart,2] <- tasksCobra$DESCRIPTION[taskStart] #Description
  tasksRaven[taskStart,3] <- tasksCobra$`SHOULD FAIL`[taskStart] #Should fail
  tasksRaven[taskStart,17] <- tasksCobra$COMMENTS[taskStart] #comments
}
#Converting metIndex
for (i in 1:nrow(metIndex)){
  metIndex$Met_ID[i] <- paste0(substr(metIndex$Met_ID[i],1,8),substr(metIndex$Met_ID[i],10,10))
}


wb = openxlsx::createWorkbook()
addWorksheet(wb,sheetName = "TASKS")
addWorksheet(wb,sheetName = "met_mapping")
writeData(wb,sheet = "TASKS",x = tasksRaven)
writeData(wb,sheet = "met_mapping",metIndex)

for (i in 1:length(taskIndEnd)){
  addStyle(wb,sheet = "TASKS",style = createStyle(border = "Bottom"),rows = taskIndEnd[i]+1,cols = 1:ncol(tasksRaven))
  
}
addStyle(wb,sheet = "TASKS",style = createStyle(border = "Bottom",borderStyle = "double"),rows = 1,cols = 1:ncol(tasksRaven))
openxlsx::saveWorkbook(wb,outputFilePath,overwrite = TRUE)



#openxlsx::::write.xlsx(tasksRaven,outputFilePath,sheetName = "TASKS",col.names = TRUE,row.names = FALSE,showNA = FALSE)
