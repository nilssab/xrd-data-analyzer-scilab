// script to convert proprietary XRD measurement data into .csv
bRes = csvDefault("eol", "linux")
files = ls(pwd()) 

for i = grep(files, '.raw')                         //set file filter           1
    [path,fname,extension]=fileparts(files(i))
    printf("\n" + fname + "\n")
    [fileD, Error] = mopen(files(i), 'rb')
    mseek(2958)
    scanSpd = mget(1, 'f', fileD)
    scanSpd_sec = scanSpd/60
    printf("Scan Speed = %f \n", scanSpd)
    startAngle = mget(1, 'f', fileD)
    printf("Start 2theta = %f \n", startAngle)
    stopAngle = mget(1, 'f', fileD)
    printf("Stop 2theta = %f \n", stopAngle)
    angleStep = mget(1, 'f', fileD)
    angleStep = angleStep*1000
    angleStep = round(angleStep)
    angleStep = angleStep/1000
    printf("2theta angle step = %f \n", angleStep)
    
    twoTheta = startAngle : angleStep : stopAngle
    timePerStep = angleStep/scanSpd_sec
    mseek(3158)
    Counts = mget(length(twoTheta), 'f', fileD)./timePerStep
    Data = [twoTheta', Counts']
    mclose(fileD)
    
    csvWrite(Data, "../"+ fname + ".csv", ',', '.',5) //set output path         2
end
