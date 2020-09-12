function [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)

    // Data:                a 2*x matrix with raw XRD data
    // d_tsh_modifier:      +-diff threshold to count as peak (change to count more/less peaks)
    peak_skirt = floor(1/((Data(1,2)-Data(1,1))*s_peakskirt_mod)) //peak skirt halfwidth (change to match y0 better)
    //FFTcutOFF:            FFT cut off range
    //plot_ok:              0 = no plots or message boxes, 1 = plots and messageboxes
    
    V_x_c = 0
    V_width = 0 
    V_Area = 0 
    V_y = 0
    y0_mean = 0
    
    Kalpha2 = 0.179285//Cu wavelengths
    Kalpha1 = 0.1788965
    Kalpha = (Kalpha2+2*Kalpha1)/3
    
    //FFT filter data
    FFTdata = fft(Data(2,:))
    a = real(FFTdata)
    b = imag(FFTdata)
    a(find(abs(a)<mean(abs(a))*FFTcutOFF))=0
    b(find(abs(b)<mean(abs(b))*FFTcutOFF))=0
    c = complex(a, b)
    //inverse FFT to get filtered data
    Filtered = fft(c, 1)
    
    //get d and dd of Filtered data
    dFiltered = diff(Filtered)
    dFiltered = [dFiltered(1), dFiltered]
    ddFiltered = diff(dFiltered)
    ddFiltered = [ddFiltered(1), ddFiltered]
    
    if plot_ok then
        figure(0)
        clf();
        plot(Data(1,:),Data(2,:))
        plot2d(Data(1,:), Filtered)
    end
    //find all starts and ends
    starts = find(dFiltered>(d_tsh_modifier*sqrt(2*max(dFiltered))))
    ends = find(dFiltered<-(d_tsh_modifier*sqrt(2*max(dFiltered))))
    //reset variables
    peak = []
    peakstart = []
    peakend = []
    nonpeak = []
    n = 0
    done = 0
    if plot_ok then
        figure(0)
    end
    //find peaks between starts and ends
    
    while (done == 0) & length(starts) & length(ends)
        n = n+1
        peakstart(n) = min(starts) //peaks start from the first +diff above threshold
        ends = ends(find(ends>min(starts))) 
        starts = starts(find(starts>min(ends))) //+diff set for peak n is removed
        
        if length(starts) == 0 then
            done = 1    //done when all +diff sets have been removed
            peakend(n) = max(ends) //when done find last peak end
        else
            peakend(n) = max(ends(find(ends<min(starts)))) //if more peaks exist, peak end is the last -diff
                                                           //below threshold and below next peak start
                                                           
            ends = ends(find(ends>min(starts))) //-diff set for peak n is removed
        end
        peak(n) = find(ddFiltered==min(ddFiltered(peakstart(n):peakend(n))))//peak position is the lowest -ddiff 
                                                                            //in respective peak area
        
        //map the region without peaks
        //peak skirt increases range that belongs to peaks                                                                    
        if n == 1 then 
            nonpeak = 1:(peakstart(n)-peak_skirt)
        end
        if n > 1 then
            nonpeak = [nonpeak, (peakend(n-1)+peak_skirt):(peakstart(n)-peak_skirt)]
        end
        if done ==1 then
            nonpeak = [nonpeak, (peakend(n)+peak_skirt):length(Data(1,:))]
        end
        
        //plot peak positions (Data and filtered already plotted)
        if plot_ok then
            plot2d([Data(1,peak(n)),Data(1,peak(n))],[0,max(Data(2,:))])
            plot2d(Data(1,peakstart(n)),Data(2,peak(n))/8+min(Data(2,:)),-1)
            plot2d(Data(1,peakend(n)),Data(2,peak(n))/8+min(Data(2,:)),-1)
        end
    end
    
    //Promt if nr of peaks are correct------------------------------------------ 
    if plot_ok then
    peaks_ok = messagebox("Peaks OK?", "modal", "question", ["Yes" "too few" "too many" "FFT wavy" "FFT noisy" "manual"])
        if peaks_ok ~= 1 then
            if peaks_ok == 0 then
                abort
            end
            if peaks_ok == 2 then
                d_tsh_modifier = d_tsh_modifier*0.9
                [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
                [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
            end
            if peaks_ok == 3 then
                d_tsh_modifier = d_tsh_modifier*1.1
                [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
                [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
            end
            if peaks_ok == 4 then
                FFTcutOFF = FFTcutOFF*0.9
                [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
                [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
            end
            if peaks_ok == 5 then
                FFTcutOFF = FFTcutOFF*1.1
                [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
                [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
            end
            if peaks_ok == 6 then
                figure(0)
                clf();
                plot(Data(1,:),Data(2,:))
                plot2d(Data(1,:), Filtered)
                title("Click peaks to mark them. Click right of Figure: OK, left: Cancel")
                
                show_window(); //put the window on the top
                
                [b,xc,yc]=xclick();
                n = 0
                peak = []
                peakstart = []
                peakend = []
                while (min(Data(1,:)) < xc) & (xc < max(Data(1,:)))//get clicks until click outside graph
                    n = n+1
                    peak(n) = floor((xc-min(Data(1,:)))/(Data(1,2)-Data(1,1)))
                    peakstart(n) = peak(n) - floor(0.5/(Data(1,2)-Data(1,1)))
                    peakend(n) = peak(n) + floor(0.5/(Data(1,2)-Data(1,1)))
                    if peakstart(n) < 1 then
                        peakstart(n) = 1
                    end
                    if peakend(n) > length(Data(1,:)) then
                        peakend(n) = length(Data(1,:))
                    end
                    
                    plot2d([Data(1,peak(n)),Data(1,peak(n))],[0,max(Data(2,:))])
                    plot2d(Data(1,peakstart(n)),Data(2,peak(n))/8+min(Data(2,:)),-1)
                    plot2d(Data(1,peakend(n)),Data(2,peak(n))/8+min(Data(2,:)),-1)
                    [b,xc,yc]=xclick();
                end
                if xc < min(Data(1,:)) then //restart if click left of figure
                    [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
                end
                if xc > max(Data(1,:)) then //peaks OK if click right of figure
                    peaks_ok = 1
                    for i = 1:length(peak)
                        if i == 1 then 
                            nonpeak = 1:(peakstart(i)-peak_skirt)
                        end
                        if i > 1 then
                            nonpeak = [nonpeak, (peakend(i-1)+peak_skirt):(peakstart(i)-peak_skirt)]
                        end
                        if i == length(peak) then
                            nonpeak = [nonpeak, (peakend(i)+peak_skirt):length(Data(1,:))]
                        end
                    end
                end
                
            end
            if peaks_ok ~= 1 then
                [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c,V_width,V_Area,V_y, y0_mean)
            end
        end
    end
    
    //mark Kalpha 2 peaks
    Kalpha2peaks = []
    n = 1
    if length(peak) > 1 then
        for i = 2:length(peak)
            Kalpha1peak = 360.*asin(Kalpha1./(2.*(Kalpha2./(2.*sin(%pi*Data(1,floor(abs(peak(i))))./(360))))))./(%pi)
            if abs(Data(1,peak(i-1))-Kalpha1peak) < 0.05 then
                Kalpha2peaks(n) = i;
                n = n+1;
            end
        end
    end
    
    //Pearson VII prelimerary parameter evaluation------------------------------
    x = [Data(1, nonpeak)', ones(length(Data(1, nonpeak)'),1)]
    y = Data(2, nonpeak)'
    z = x\y         //find y0 by least squares on nonpeak regions 
    y0_a = z(1)
    y0_b = z(2)
    x = Data(1,:)
    y0 = x.*y0_a+y0_b
    y0_1 = min(Data(1,:))*y0_a+y0_b
    y0_2 = max(Data(1,:))*y0_a+y0_b
    
    //reset values
    x_c = []
    Area = []
    width = []
    y = y0
    n = 1
    m = 0
    for i = 1:length(peak) //find x_c, width, Area for non Kalpha2 peaks
        if i ~= Kalpha2peaks(n) then
            m = m+1
            if length(peak)-i > 0 then
                if i+1 == Kalpha2peaks(n) then  //fit using Kaplha1 peak position if next peak is Kalpha2
                    x_c(m) = Data(1,peak(i))
                else                            //convert to Kaplha peak position if peak contains Kalpha 2
                    x_c(m) = 360.*asin(Kalpha1./(2.*(Kalpha./(2.*sin(%pi*Data(1,peak(i))./(360))))))./(%pi) 
                end
            else                                //convert to Kaplha peak position if peak contains Kalpha 2
                x_c(m) = 360.*asin(Kalpha1./(2.*(Kalpha./(2.*sin(%pi*Data(1,peak(i))./(360))))))./(%pi) //Kaplha 1 peak position
            end
            
            Area(m) = inttrap(Data(1,peakstart(i):peakend(i)), Filtered(peakstart(i):peakend(i)) - y0(peakstart(i):peakend(i)))/1.35//Peak area estimate
            width(m) = (Data(1,peakend(i))-Data(1,peakstart(i)))/3                              //Peak width estimate
        else
            n = n +1
        end
    end
    
    //evaluate y with all x_c peaks by Pearson VII
    y = y_func(x, x_c, width, Area, y0_a, y0_b)
    
    if plot_ok then //plot fitted curve over XRD data
        figure(0)
        clf();
        plot(x, Data(2,:))
        plot2d(x, y)
    end
    
    //Peak Fitting--------------------------------------------------------------
    printf("Error = %f \n", inttrap(x,(Data(2,:)-y).^2)) 
    //continue_fitting = messagebox("Continue fitting?", "modal", "question", ["Yes" "No" "fit y0"])
    
    gamma_x = [] //set fit step sizes (gamma) 
    gamma_w = []
    gamma_A = []
    for i = 1:length(x_c)
        gamma_x(i) = 0.0000004*g_factor
        gamma_w(i) = width(i)/800000*g_factor
        gamma_A(i) = Area(i)/600000*g_factor
    end
    gamma_y0_1 = y0_1/300*g_factor
    gamma_y0_2 = y0_2/300*g_factor
    
    ErrorNew = inttrap(x,(Data(2,:)-y).^2)
    ErrorOld = ErrorNew*2
    iteration = 1
    done = 0
    //Main part of fitting 
    //ErrorOld/ErrorNew > 1.000001) | (ErrorOld/ErrorNew < 1
    while (~done)
        iteration = iteration + 1
        ErrorOld = ErrorNew
        for i = 1:length(x_c)//get the dE gradient for all parameters
            
            x_temp = x_c
            x_temp(i) = x_temp(i) + gamma_x(i)
            dEdx(i) = (inttrap(x,(Data(2,:)-y_func(x, x_temp, width, Area, y0_a, y0_b)).^2) - ErrorNew)
            x_temp = x_c
            x_temp(i) = x_c(i) - (gamma_x(i) * dEdx(i))
            y_temp = y_func(x, x_temp, width, Area, y0_a, y0_b)
            ErrorNew_temp = inttrap(x,(Data(2,:)-y_temp).^2)
            if  ErrorNew_temp < ErrorNew then
                x_c = x_temp
                ErrorNew = ErrorNew_temp
                gamma_x(i) = gamma_x(i) * 5
            else
                gamma_x(i) = gamma_x(i) * 0.5
            end
            
            w_temp = width
            w_temp(i) = w_temp(i) + gamma_w(i)
            dEdw(i) = (inttrap(x,(Data(2,:)-y_func(x, x_c, w_temp, Area, y0_a, y0_b)).^2) - ErrorNew)
            w_temp = width
            w_temp(i) = width(i) - (gamma_w(i) * dEdw(i))
            y_temp = y_func(x, x_c, w_temp, Area, y0_a, y0_b)
            ErrorNew_temp = inttrap(x,(Data(2,:)-y_temp).^2)
            if  ErrorNew_temp < ErrorNew then
                width = w_temp
                ErrorNew = ErrorNew_temp
                gamma_w(i) = gamma_w(i) * 5
            else
                gamma_w(i) = gamma_w(i) * 0.5
            end
            
            A_temp = Area
            A_temp(i) = A_temp(i) + gamma_A(i)
            dEdA(i) = (inttrap(x,(Data(2,:)-y_func(x, x_c, width, A_temp, y0_a, y0_b)).^2) - ErrorNew)
            A_temp = Area
            A_temp(i) = Area(i) - (gamma_A(i) * dEdA(i))
            y_temp = y_func(x, x_c, width, A_temp, y0_a, y0_b)
            ErrorNew_temp = inttrap(x,(Data(2,:)-y_temp).^2)
            if  ErrorNew_temp < ErrorNew then
                Area = A_temp
                ErrorNew = ErrorNew_temp
                gamma_A(i) = gamma_A(i) * 5
            else
                gamma_A(i) = gamma_A(i) * 0.5
            end
            
        end
        
        y0_1_temp = y0_1
        y0_1_temp = y0_1_temp + gamma_y0_1
        x_y0 = [min(Data(1,:)), max(Data(1,:))]
        y_y0 = [y0_1_temp, y0_2]
        z_y0 = x_y0\y_y0         //find y0 by least squares on nonpeak regions 
        y0_a_temp = z_y0(1)
        y0_b_temp = z_y0(2)
        dEdy0_1 = (inttrap(x,(Filtered - y_func(x, x_c, width, Area, y0_a_temp, y0_b_temp)).^2) - ErrorNew)
        y0_1_temp = y0_1 - (gamma_y0_1 * dEdy0_1)
        x_y0 = [min(Data(1,:)), max(Data(1,:))]
        y_y0 = [y0_1_temp, y0_2]
        z_y0 = x_y0\y_y0         //find y0 by least squares on nonpeak regions 
        y0_a_temp = z_y0(1)
        y0_b_temp = z_y0(2)
        y_temp = y_func(x, x_c, width, Area, y0_a_temp, y0_b_temp)
        ErrorNew_temp = inttrap(x,(Data(2,:)-y_temp).^2)
        if  ErrorNew_temp < ErrorNew then
            y0_1 = y0_1_temp
            y0_a = y0_a_temp 
            y0_b = y0_b_temp
            ErrorNew = ErrorNew_temp
            gamma_y0_1 = gamma_y0_1 * 5
        else
            gamma_y0_1 = gamma_y0_1 * 0.5
        end
        
        y0_2_temp = y0_2
        y0_2_temp = y0_2_temp + gamma_y0_2
        x_y0 = [min(Data(1,:)), max(Data(1,:))]
        y_y0 = [y0_1, y0_2_temp]
        z_y0 = x_y0\y_y0         //find y0 by least squares on nonpeak regions 
        y0_a_temp = z_y0(1)
        y0_b_temp = z_y0(2)
        dEdy0_2 = (inttrap(x,(Filtered - y_func(x, x_c, width, Area, y0_a_temp, y0_b_temp)).^2) - ErrorNew)
        y0_2_temp = y0_2 - (gamma_y0_2 * dEdy0_2)
        x_y0 = [min(Data(1,:)), max(Data(1,:))]
        y_y0 = [y0_1, y0_2_temp]
        z_y0 = x_y0\y_y0         //find y0 by least squares on nonpeak regions 
        y0_a_temp = z_y0(1)
        y0_b_temp = z_y0(2)
        y_temp = y_func(x, x_c, width, Area, y0_a_temp, y0_b_temp)
        ErrorNew_temp = inttrap(x,(Data(2,:)-y_temp).^2)
        if  ErrorNew_temp < ErrorNew then
            y0_2 = y0_2_temp
            y0_a = y0_a_temp 
            y0_b = y0_b_temp
            ErrorNew = ErrorNew_temp
            gamma_y0_2 = gamma_y0_2 * 5
        else
            gamma_y0_2 = gamma_y0_2 * 0.5
        end

        y = y_func(x, x_c, width, Area, y0_a, y0_b)//evaluate new y
        
        ErrorNew = inttrap(x,(Data(2,:)-y).^2)
        
        printf("Error = %f ", ErrorNew)
        printf("peak = %f ", x_c)
        printf("\n") 
        //continue_fitting = messagebox("Continue fitting?", "modal", "question", ["Yes" "No" "fit y0"])
        
        
        if iteration >= iterations then
            done = 1
            if plot_ok then //plot new fitted curve over XRD data
                figure(0)
                clf()
                plot(x, Data(2,:))
                plot2d(x, y)
            end
            if plot_ok then //ask if fitting is satisfactory                    1       2       3       4               5           6
                fitting_ok = messagebox("Fitting OK?", "modal", "question", ["Yes" "gamma+" "gamma-" "manual y0" "use filtered" "Abort"])
            end
            if fitting_ok == 4 then
                done = 0
                iteration = 1
                
                figure(0)
                
                title("Click two y0 points to manually fit to")
                
                show_window(); //put the window on the top
                
                [b,xc,yc]=xclick();
                y0_1 = yc
                [b,xc,yc]=xclick();
                y0_2 = yc
                x_y0 = [min(Data(1,:)), 1; max(Data(1,:)),1]
                y_y0 = [y0_1; y0_2]
                z_y0 = x_y0\y_y0         //find y0 by least squares on nonpeak regions 
                y0_a = z_y0(1)
                y0_b = z_y0(2)
                printf("y0_ 1 2 = %f, %f. z = %f, %f.\n",y0_1, y0_2, z_y0(1), z_y0(2))
                y = y_func(x, x_c, width, Area, y0_a, y0_b)//evaluate new y
                if plot_ok then //plot new fitted curve over XRD data
                    figure(0)
                    clf()
                    plot(x, Data(2,:))
                    plot2d(x, y)
                end
            end
        end
        
        
    end
    //error is evaluated by inttrap(x,(Data(2,:)-y).^2)

    
    
    if plot_ok then 
        if fitting_ok == 6 then
            [V_x_c, V_width, V_Area, V_y, y0_mean] = return(0,0,0,0)
        end
        if fitting_ok == 0 then
            abort
        end
        if fitting_ok == 2 then
            g_factor = g_factor*3
            [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
            [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
        end
        if fitting_ok == 3 then
            g_factor = g_factor*0.2
            [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
            [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
        end
        
        if fitting_ok == 5 then
            Data = [x; Filtered]
            [V_x_c, V_width, V_Area, V_y, y0_mean] = find_peaks(Data, d_tsh_modifier, s_peakskirt_mod, FFTcutOFF, plot_ok, g_factor, iterations)
            [V_x_c, V_width, V_Area, V_y, y0_mean] = return(V_x_c, V_width, V_Area, V_y, y0_mean)
        end
    end
    
    V_y = y //return values and end
    V_x_c = x_c
    V_width = width
    V_Area = Area
    y0_mean = (y0_1+y0_2)/2
endfunction
