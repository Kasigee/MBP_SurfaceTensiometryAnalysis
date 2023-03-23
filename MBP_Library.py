import numpy as np
import pandas as pd
from scipy import stats
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
from adjustText import adjust_text
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams, cycler
from matplotlib.patches import Patch
from scipy.fft import rfft,rfftfreq
plt.rc('xtick', labelsize=14)     
plt.rc('ytick', labelsize=12)

# %matplotlib inline #not used now this has been turned into a library.

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 10        }

mpl.rcParams['figure.dpi'] = 200

plt.rc('font', **font)

def removesuffix(self: str, suffix: str, /) -> str:
    # suffix='' should not call self[:-0].
    if suffix and self.endswith(suffix):
        return self[:-len(suffix)]
    else:
        return self[:]
    

# class MyListState:
#     def __init__(self):
#         self.SectionTimeStart = []
#         self.CUMULATIVE_RollingMAX = []
#         self.CUMULATIVE_AvgMAX = []
#         self.CUMULATIVE_AvgPeak = []
#         self.CUMULATIVE_RollingMAX_err = []
#         self.CUMULATIVE_AvgPeak_err = []

# def Initialise_variables(Baseline,Sensitivity,EtOHST,solvent_density):
#     Baseline=Baseline # Voltage of Pressure sensor in air. Note baseline weeks later gives: Average Peak: 1.86373232589543 +/- 0.000393771406842909 (STDEV) or 9.722750786244667e-06 (STERR) or 0.021128109512921687 % (STDEV%) OR 0.0005216817163684369 % (STERR%) -- Filename: BaselineToday_20221129_120031
#     Sensitivity=Sensitivity #Pa/V  --> This is for Pressure sensor "2" - I'll find the actual code on the side to identify it properly. Connected currently to the low dead volume block. Calibrated via hydrostatic pressure measurements in triplicate.
#     EtOHST=EtOHST #N/m
#     solvent_density=solvent_density #kg.m^-3 EtOH
#     return Baseline, Sensitivity, EtOHST, solvent_density

SectionTimeStart = []
CUMULATIVE_RollingMAX = []
CUMULATIVE_AvgMAX = []
CUMULATIVE_AvgPeak = []
CUMULATIVE_RollingMAX_err = []
CUMULATIVE_AvgMAX_err = []
CUMULATIVE_AvgPeak_err = []
# CUMULATIVE_PackedMAX=[] #This is not used yet, but is kept for future reference.

def refresh_states():
    # SectionTimeStart = []
    # CUMULATIVE_RollingMAX = []
    # CUMULATIVE_AvgMAX = []
    # CUMULATIVE_AvgPeak = []
    # CUMULATIVE_RollingMAX_err = []
    # CUMULATIVE_AvgPeak_err = []
    global SectionTimeStart, CUMULATIVE_RollingMAX, CUMULATIVE_AvgMAX, CUMULATIVE_AvgPeak, CUMULATIVE_RollingMAX_err, CUMULATIVE_AvgMAX_err, CUMULATIVE_AvgPeak_err


def peak_analysis(Time_between_bubbles,Freq,split, y, i, length):
    AVERAGEPEAKarray=[]
    AVERAGEPEAKerrarray=[]
    peaks, _ = find_peaks(y, distance=Time_between_bubbles*Freq)
    peaks=peaks+((i-1)*int(length/split))
    peaks_min, _ = find_peaks(-y, distance=Time_between_bubbles*Freq)
    peaks_min=peaks_min+((i-1)*int(length/split))
    local_y_max=y[peaks].max()
    #print("local max =",local_y_max)
    for idx in (y[peaks].index):
        if y[idx] < (local_y_max*0.995):
            peaks = np.delete(peaks, np.where(peaks == idx))
    LengthPeak=len(peaks)
    AVERAGEPeak=y[peaks].mean()
    STDPeak=y[peaks].std()
    StErrPeak=STDPeak/(LengthPeak*(1/2))
    PercentPeak=(STDPeak/AVERAGEPeak)*100
    PercentPeak2=(StErrPeak/AVERAGEPeak)*100
    print("Average Peak:",AVERAGEPeak,"+/-",STDPeak,"or",StErrPeak,"or",PercentPeak,"% OR", PercentPeak2,"%")
    AVERAGEPEAKarray.append(AVERAGEPeak)
    AVERAGEPEAKerrarray.append(STDPeak)
    LengthPeakmin=len(peaks_min)
    # print(LengthPeakmin)
    AVERAGEPeakmin=Baseline
    STDPeakmin=y[peaks_min].std()
    StErrPeakmin=STDPeakmin/(LengthPeakmin*(1/2))
    PercentPeakmin=(STDPeakmin/AVERAGEPeakmin)*100
    PercentPeak2min=(StErrPeakmin/AVERAGEPeakmin)*100
    print("Average Peak min:",AVERAGEPeakmin,"+/-",StErrPeakmin,"or",PercentPeakmin,"% OR", PercentPeak2min,"%")
    P4=AVERAGEPeak*Sensitivity
    E4=StErrPeak*Sensitivity
    return AVERAGEPeak, STDPeak, StErrPeak, PercentPeak, PercentPeak2, AVERAGEPeakmin, STDPeakmin, StErrPeakmin, PercentPeakmin, PercentPeak2min, P4, E4

def LargestSmallestPoint(y, Sensitivity, EtOHST, solvent_density, Baseline):
    Largest=y.nlargest(n=20) #Do this via peaks instead? This checks the largest value 20 points.
    LargestSTD=Largest.std()
    Largest=Largest.mean() #This is the average of the 20 largest points.
    Smallest=y.nsmallest(n=20)
    Smallest=Smallest.mean()
    Smallest=Baseline
    V_Diff=Largest-Smallest
    # P1=Sensitivity*Largest
    # P_Diff=Sensitivity*V_Diff
    CalculateHydrostaticPressure(solvent_density, 0.001) #Assuming a 0.001m depth of the needle in the solvent
    print("V from the average of the 20 largest points is calculated to be", '\033[1m', Largest,'\033[0m',"+/-",LargestSTD,"V. Hmm, it takes the y-value from the rolling average, not the actual data. This solves the problem of not gathering all the points from the same peak, that could happen.")
    return Largest, LargestSTD

# def CalculateHydrostaticPressure(solvent_density, Depth, altitude):
#     radiusOfEarth=6371009 #m
#     graviational_acceleration=9.80665*(radiusOfEarth/(radiusOfEarth+altitude))**2 #Assuming on earth with a gravitational acceleration of 9.80665m/s^2 at sea level.
#     print("Graviational acceleration is", graviational_acceleration, "m/s^2 (assuming an altitude of", altitude, "m above sea level)")
#     HP=solvent_density*Depth*graviational_acceleration
#     print("Hydrostatic pressure is:",HP,"Pa (assuming a ", Depth, "m depth of the needle in the solvent with density", solvent_density, "kg/m^3)")
#     return HP

#other changes are going to have a greater affect then the change in gravitational constant.
def CalculateHydrostaticPressure(solvent_density, Depth):
    graviational_acceleration=9.80665 #Assuming on earth with a gravitational acceleration of 9.80665m/s^2 at sea level.
    # print("Graviational acceleration is", graviational_acceleration, "m/s^2")
    HP=solvent_density*Depth*graviational_acceleration
    print("Hydrostatic pressure is:",HP,"Pa (assuming a ", Depth, "m depth of the needle in the solvent with density", solvent_density, "kg/m^3)")
    return HP

def convert_STtoP(CapillaryRadius, ST):
    P_Diff=(ST*2)/CapillaryRadius
    print("Expected pressure from the surface tension (ignoring hydrostatic pressure) is", P_Diff, "Pa")
    HP=CalculateHydrostaticPressure(997, 0.01) #Assuming a 0.01m depth of the needle in water . Density of water is 997kg/m^3.
    P_error_ignoring_hydrostatic_pressure=HP/(HP+P_Diff)
    print("At a depth of 1cm in water, the hydostatic pressure is", HP, "Pa", "so the total pressure is", P_Diff+HP, "Pa. The error in the pressure measurement is", P_error_ignoring_hydrostatic_pressure*100, "%")
    return P_Diff

def convert_VtoP(Sensitivity, V_Diff):
    P_Diff=Sensitivity*V_Diff
    print("P is", P_Diff, "Pa")
    return P_Diff

def convert_PtoST(CapillaryRadius, P_Diff, solvent_density, Depth):
    HP=CalculateHydrostaticPressure(997, 0.01) #Assuming a 0.01m depth of the needle in water. Density of water is 997kg/m^3.
    ST=CapillaryRadius/((P_Diff-HP)*2)
    print("ST is", ST, "mN/m")
    return ST

def CalibrateRadius(ST,P_Diff): #ST is in mN/m and P_Diff is in Pa
    HP=CalculateHydrostaticPressure(997, 0.01) #Assuming a 0.01m depth of the needle in water. Density of water is 997kg/m^3.
    CapillaryRadius=2*ST/(P_Diff-HP)
    return CapillaryRadius

def CalculateContactAngle(k,DP1,DP2,ST):
    contactAngle=arccos((k*DP1*DP2)/(2*ST(DP2-DP1)))
    print("Contact angle is", contactAngle)
    return contactAngle

def PackingData(i,length,split,packet_size,RawY,RawX,averages,time):
    # averages = []
    # time = []
    for j in range(((i-1)*int(length/split)), i*int(length/split), packet_size):
        # Get the current packet of data
        packet = RawY[j:j+packet_size]
        packetTime = RawX[j:j+packet_size]
        # Calculate the average of the packet
        average = sum(packet) / len(packet)
        averageTime = (sum(packetTime) / len(packetTime))/(60)
        # Add the average to the list
        averages.append(average)
        time.append(averageTime)
    return averages, time

def RollingAvgData(Data,Freq,avgPerSec,split,FileName, i, length,time,averages):
    RollingSTD = Data.rolling(Freq).std() ## This is the standard deviation of the rolling average. It might need to have its interval changed to Freq/avgPerSec???
    Data1= Data
    Data = Data.rolling(int(Freq/avgPerSec)).mean() # Maybe take legit max/min and minus the noise.
    print("If using times 3, need to ensure there is longer time left in the minimum state, otherwise it fails.")
    x=Data[Data.columns[0]].astype(float)
    y=Data[Data.columns[1]].astype(float)
    print("length of x is = ", len(x), " and length of Rolling y is = ",len(y))
    ySTD=RollingSTD[RollingSTD.columns[1]].astype(float)
    ymin=y-ySTD
    ymax=y+ySTD
    x1=Data1[Data1.columns[0]].astype(float)
    y1=Data1[Data1.columns[1]].astype(float)
    plt.plot(x/60,y,color='red')
    plt.plot(time,averages,color='black', linestyle='--')
    return x, y, ySTD, ymin, ymax, x1, y1
    
def rollingMAXMIN_analysis(Time_between_bubbles,SwitchFreq,Freq,split,FileName, y, i, length):
    SwitchFreq=int((Time_between_bubbles*1.5)*Freq)
    RollingMax= y.rolling(SwitchFreq).max().rolling(20).mean() #FIX THIS; it does seem to be an improvement though.
    RollingMaxSTD= y.rolling(SwitchFreq).max().rolling(20).std() #FIX THIS; it does seem to be an improvement though.
    RollingMin= y.rolling(SwitchFreq).min().rolling(20).mean()
    RollingDiff=(RollingMax-RollingMin).mean()
    RollingDiff=(RollingMax-Baseline).mean()
    RollingSTD=(RollingMax-RollingMin).std()
    PercentErr=RollingSTD/RollingDiff
    P2=Sensitivity*RollingMax
    E2=Sensitivity*RollingMaxSTD
    V2=RollingMax.mean()
    E2=RollingMax.std()
    RollingMaxAvg=RollingMax.mean()
    RollingMinAvg=RollingMin.mean()
    RollingMaxStd=RollingMax.std()
    RollingMinStd=RollingMin.std()
    RollingAVGDiff=RollingMaxAvg-Baseline
    RollingAVGStd=RollingMaxStd+RollingMaxStd
    PercentErr2=RollingAVGStd/RollingAVGDiff
    P3=Sensitivity*RollingMaxAvg
    E3=Sensitivity*RollingMaxStd
    return RollingMax, RollingMaxSTD, RollingMin, RollingDiff, RollingSTD, PercentErr, P2, E2, V2, E2, RollingMaxAvg, RollingMinAvg, RollingMaxStd, RollingMinStd, RollingAVGDiff, RollingAVGStd, PercentErr2, P3, E3

def plot_rolling(x,RollingMax,RollingMin):
    plt.plot(x/60,RollingMax)
    plt.plot(x/60,RollingMin)
    plt.ylabel('ΔV(V)')
    plt.xlabel('Time(min)')

def FindTotalTime(Data,frequency):
    length=len(Data)
    TotalTime=length/frequency
    print("Total time is", TotalTime, "seconds")
    return TotalTime

def DataOverTime(Time,RollingMax2,Largest,AVERAGEPeak,RollingMaxSTD2,LargestSTD,STDPeak):
    SectionTimeStart.append(Time)
    CUMULATIVE_RollingMAX.append(RollingMax2)    
    CUMULATIVE_AvgMAX.append(Largest)
    CUMULATIVE_AvgPeak.append(AVERAGEPeak)
    CUMULATIVE_RollingMAX_err.append(RollingMaxSTD2)    
    CUMULATIVE_AvgMAX_err.append(LargestSTD)
    CUMULATIVE_AvgPeak_err.append(STDPeak)
    # CUMULATIVE_PackedMAX.append(PackedMax) #This is not used yet, but is kept for future reference.
    print("Cumulative SectionTimeStart Series is", SectionTimeStart)
    print("Cumulative RollingMax is", CUMULATIVE_RollingMAX)
    print("Cumulative AvgMax is", CUMULATIVE_AvgMAX)
    print("Cumulative AvgPeak is", CUMULATIVE_AvgPeak)
    print("Cumulative RollingMaxSTD is", CUMULATIVE_RollingMAX_err)
    print("Cumulative LargestSTD is", CUMULATIVE_AvgMAX_err)
    print("Cumulative AvgPeakSTD is", CUMULATIVE_AvgPeak_err)
    return SectionTimeStart, CUMULATIVE_RollingMAX, CUMULATIVE_AvgMAX, CUMULATIVE_AvgPeak, CUMULATIVE_RollingMAX_err, CUMULATIVE_AvgMAX_err, CUMULATIVE_AvgPeak_err


#If you see blatently wrong maths compared to the variable, I was likely testing something and forgot to change it back. If you see subtlely wrong maths, I was likely just wrong.    
def plot_pressure(FileName,Freq,split,Time_between_bubbles):
    # SectionTimeStart=[]
    # CUMULATIVE_RollingMAX=[]
    # CUMULATIVE_AvgMAX=[]
    # CUMULATIVE_AvgPeak=[]
    # CUMULATIVE_RollingMAX_err=[]
    # CUMULATIVE_AvgPeak_err=[]
    # CUMULATIVE_PackedMAX=[] #This is not used yet, but is kept for future reference.
    # return SectionTimeStart, CUMULATIVE_RollingMAX, CUMULATIVE_AvgMAX, CUMULATIVE_AvgPeak, CUMULATIVE_RollingMAX_err, CUMULATIVE_AvgPeak_err

    # SectionTimeStart=[]
    avgPerSec=10
    packet_size=int(Freq/avgPerSec)
    SwitchFreq=int((Time_between_bubbles*1.5)*Freq)
    Output=pd.DataFrame()
    
    Data_file = pd.read_csv(r'{}.csv'.format(FileName), on_bad_lines='skip', delimiter=',')
    Data_file=Data_file[4:] #Removes header info
    Data_file=Data_file[5*Freq:] #Removes first 5 seconds of data
    
    for i in range(1,split+1):
        Data = Data_file.reset_index(drop=False) #(Re)initialises the Data.
        Data.drop(columns=Data.columns[-1], axis=1, inplace=True)
        RawX=Data[Data.columns[0]].astype(float)
        RawY=Data[Data.columns[1]].astype(float)
        length=len(RawY) 
        
        # Initialize a list to store the packet averages
        # Iterate over the data in packets of size `packet_size`
        print("Packing data into packets of size", packet_size, "for", i, "of", split, "packets")
        averages = []
        time = []
        averages,time=PackingData(i,length,split,packet_size,RawY,RawX,averages,time)
        # print("average array:", averages)
        # print("time array:", time)
        

        length=len(Data)
        Data=Data[((i-1)*int(length/split)):i*int(length/split)]
        print("Rolling data")
        x, y, ySTD, ymin, ymax, x1, y1 = RollingAvgData(Data,Freq,avgPerSec,split,FileName, i, length,time,averages)

        plt.legend(['Rolling','Packed'])
        plt.text(0.95, 0.95, 'Sensitivity: {}\nEtOHST: {}\nsolvent_density: {}\nFreq: {}\nTime_between_bubbles: {}\nBaseline: {}\nsplit: {}\navgPerSec: {}'.format(Sensitivity, EtOHST, solvent_density, Freq, Time_between_bubbles, Baseline, i, avgPerSec), ha='right', va='top', transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.5))

        print("Peak analysis")
        AVERAGEPeak, STDPeak, StErrPeak, PercentPeak, PercentPeak2, AVERAGEPeakmin, STDPeakmin, StErrPeakmin, PercentPeakmin, PercentPeak2min, P4, E4 = peak_analysis(Time_between_bubbles,Freq,split,y,i,length)

        print("Largest and smallest test.")
        Largest, LargestSTD = LargestSmallestPoint(y, Sensitivity, EtOHST, solvent_density, Baseline)
        print("Rolling max and min analysis")
        RollingMax, RollingMaxSTD, RollingMin, RollingDiff, RollingSTD, PercentErr, P2, E2, V2, E2, RollingMaxAvg, RollingMinAvg, RollingMaxStd, RollingMinStd, RollingAVGDiff, RollingAVGStd, PercentErr2, P3, E3 = rollingMAXMIN_analysis(Time_between_bubbles,SwitchFreq,Freq,split,FileName, y, i, length)
        plot_rolling(x,RollingMax,RollingMin)

        RollingMax2=RollingMax.mean()
        RollingMaxSTD2=RollingSTD.mean()
        print("V rolling avg is calculated to be", '\033[1m', RollingMax2,'\033[0m',"+/-",RollingMaxSTD2,"V")
        #print("V avg of rolling max  is calculated to be", '\033[1m', RollingMaxAvg,'\033[0m',"+/-",RollingMaxStd,"V")
       
        plt.show()
        SectionTimeStart, CUMULATIVE_RollingMAX, CUMULATIVE_AvgMAX, CUMULATIVE_AvgPeak, CUMULATIVE_RollingMAX_err, CUMULATIVE_AvgMAX_err,  CUMULATIVE_AvgPeak_err = DataOverTime(time[0],RollingMax2,Largest,AVERAGEPeak,RollingMaxSTD2,LargestSTD,STDPeak)


    ## Saves data to a csv

        # Convert averages and time lists to Pandas Series
        averages_series = pd.Series(averages)
        time_series = pd.Series(time)
        # Concatenate the dataframes vertically
        x_processed=x.dropna()
        y_processed=y.dropna()
        time_series_processed=time_series.dropna()
        averages_series_processed=averages_series.dropna()

        Output0 = pd.concat( [x_processed,y_processed,time_series_processed,averages_series_processed], axis=1)
        Output = Output.append(Output0)
    print("Converting to csv")
    Output2 = pd.DataFrame({'SectionTimeStart':SectionTimeStart,'MaxFromRollingAvg': CUMULATIVE_RollingMAX, 'LargestPressureInSplit': CUMULATIVE_AvgMAX, 'AveragePeak(rolling)': CUMULATIVE_AvgPeak, 'MaxFromRollingAvg_err': CUMULATIVE_RollingMAX_err, 'AverageMax_err':CUMULATIVE_AvgMAX_err,'AveragePeak(rolling)_err': CUMULATIVE_AvgPeak_err})
    Output.columns = ['Time(rolling)', 'Voltage(rolling)','Time(packed)','Voltage(packed)']
    Output2.columns = ['SectionTimeStart(min)','MaxFromRollingAvg', 'LargestPressureInSplit','AveragePeak(rolling)', 'MaxFromRollingAvg_err', 'AverageMax_err','AveragePeak(rolling)_err']
    
    # Save the resulting dataframe to a CSV file
    Output.to_csv('{}_processed.csv'.format(FileName), index=False)
    Output2.to_csv('{}_DataOverTime_{}.csv'.format(FileName,split), index=False)
    x = Output2['SectionTimeStart(min)']
    y1 = Output2['MaxFromRollingAvg']
    y2 = Output2['LargestPressureInSplit']
    y3 = Output2['AveragePeak(rolling)']
    MaxFromRollingAvg_err = Output2['MaxFromRollingAvg_err']
    AverageMax_err = Output2['AverageMax_err']
    AveragePeak_rolling_err = Output2['AveragePeak(rolling)_err']
    # ax = Output2.plot(x='SectionTimeStart(min)', y=['MaxFromRollingAvg', 'LargestPressureInSplit','AveragePeak(rolling)'], kind="line", title="{} Data Over Time".format(FileName))
    # ax = Output2.errorbar(x='SectionTimeStart(min)', y=['MaxFromRollingAvg', 'AveragePeak(rolling)'], yerr=['MaxFromRollingAvg_err', 'AveragePeak(rolling)_err'], capsize=5, capthick=2)
    # Create a figure and axes object
    fig, ax = plt.subplots()

    # Create the scatter plot
    ax.scatter(x, y1, label='MaxFromRollingAvg')
    ax.scatter(x, y2, label='LargestPressureInSplit')
    ax.scatter(x, y3, label='AveragePeak(rolling)')

    # Add error bars
    ax.errorbar(x, y1, yerr=MaxFromRollingAvg_err, fmt='none', capsize=5, capthick=2, ecolor='black')
    ax.errorbar(x, y1, yerr=MaxFromRollingAvg_err, fmt='none', capsize=5, capthick=2, ecolor='black')
    ax.errorbar(x, y3, yerr=AveragePeak_rolling_err, fmt='none', capsize=5, capthick=2, ecolor='black')

    # Add title, x-label, y-label, and legend
    ax.set_title("{} Data Over Time".format(FileName))
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Voltage (V)")
    ax.legend()
    # ax.set_xlabel("Time (min)")
    # ax.set_ylabel("Voltage (V)")
    # ax.legend(['MaxFromRollingAvg', 'LargestPressureInSplit','AveragePeak(rolling)'])
    plt.show()
    # print(AVERAGEPEAKarray)
