//clear
//exec("y_func.sce path")
//ecec("find_peaks.sce path") //1.0002% limit is set in here

bRes = csvDefault("eol", "linux")

files = ls(pwd())
files = gsort(files, 'lr','i')
//peak search parameters
deltaTsh = 1                //raise if too many peaks are found                             1
peakskirt_mod = 7           //raise if peaks close to border, lower if y0 is high           2
FFTcutOFF = 2             //raise if noisy, lower if wavy                                 3
//fitting parameters
g_factor = 1                //lower if divergence is high, lower if fitting is insufficient 4
iterations = 40

n = 1

PeakData = []
max_peaks = 0
for i = grep(files, '.csv') //set filename filter                                           5
    [path,fname,extension]=fileparts(files(i))
    printf("\n" + fname + "\n")
    Data = csvRead(files(i))'
    if max(Data) > 10
        [peaks, widths, area, y, y0] = find_peaks(Data,deltaTsh,peakskirt_mod,FFTcutOFF,1,g_factor,iterations)
        if length(peaks) > max_peaks then
            max_peaks = length(peaks)
        end
    else
        peaks = 0
    end
    
    printf("peak = %f \n", peaks)
    printf("widths = %f \n", widths)
    peakData = fname + "," + string(y0) + ","
    for m = 1:length(peaks)
        peakData = peakData + string(peaks(m)) + ","//comment away unwanted fields          6
        peakData = peakData + string(widths(m)) + "," 
        peakData = peakData + string(area(m)) + ","
        if length(y) > 1 then
            peakData = peakData + string(y(floor((peaks(m)-min(Data(1,:)))/(Data(1,2)-Data(1,1)))))
        end
        
        if m ~= length(peaks) then
            peakData = peakData + ",,"
        end
    end
    
    mkdir("report")                             //create folder if used
    [fd, err] = mopen("report/report.csv", 'a') //set report output path                    7
    mputl(peakData, fd)
    mclose(fd)
    mkdir("fitdata")
    csvWrite([Data(1,:)',y'],"fitdata/" + fname + "_fitting.csv")
    n = n+1
end



