print("""
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
\tICoFit 2.0.1
\tUpdated: 2019/12/16 by Shun Tokuda
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
""")

""" Import packages """
print("Importing packages.")
import numpy as np
import pandas as pd
import sys
from scipy.optimize import curve_fit
import glob
print("  All packages were successfully imported.")

## Define fitting calculation

def StretchedExp_fit(filepath):  
    # Import the excel file
    raw_data = pd.read_excel(filepath)
    
    # names of each parameter corresponding to the paper
    def func_ICF(CDT, sigma2, A, tauf, taus, beta):
        CD =  sigma2 * (A*np.e**(-CDT/tauf)+(1-A)*np.e**(-(CDT/taus)**beta))**2 * (tauf<taus)
        return CD

    # Initial values of fitting parameters
    initial=[1,1,200,600,0.1] #sigma2, A,tauf[us],taus[us],beta
    print("  Initial values for the fitting")
    print("    Default:\n")
    print("\t sigma2 \t 1")
    print("\t A \t\t 1")
    print("\t tau_fast[µs] \t 200")
    print("\t tau_slow[µs] \t 600")
    print("\t beta \t\t 0.1\n")
    print("    Change from defaut? (y/n)")
    ini_change_determined = False
    while ini_change_determined == False:
        ini_change = input("      >> ")
        if ini_change == "n":
            ini_change_determined = True
        elif ini_change == "y":
            print("    Input new initial values:\n")
            initial[0] = float(input("\t sigma2 \t>> "))
            initial[1] = float(input("\t A \t\t>> "))
            initial[2] = float(input("\t tau_fast [µs] \t>> "))
            initial[3] = float(input("\t tau_slow [µs] \t>> "))
            initial[4] = float(input("\t beta \t\t>> "))
            ini_change_determined = True
        else:
            print("    Anser with \"y\" or \"n\"!")
            
    # Boundary condition
    print("\n  Boundary conditions:\n")
    print("\t 0 <\tsigma2\t< 1")
    print("\t 0 <\tA\t< 1")
    print("\t 0 <\ttau_fast")
    print("\t 0 <\ttau_slow")
    print("\t 0 <\tbeta\t< 1")
    
    # Constrain
    print("\n  Constrain:\n")
    print("\t tau_fast < tau_slow")
    
    # Fitted parameters and variances are collected in:
    popts=pd.DataFrame()
    sd=pd.DataFrame()
    r2=[]

    # delay time
    CDT = pd.DataFrame(raw_data.loc[1,'Correlation Delay Times[1] (µs)':'Correlation Delay Times[192] (µs)'])
    
    # Report the number of data to be fitted
    print("\n  ",len(raw_data),"data will be fitted.")
    
    # Fit
    for r in np.arange(0,len(raw_data)):
        print("\n  ",r+1,"/",len(raw_data),end="")
        try:
            # Make DataFrame containning CDT and CD
            CD = pd.DataFrame(raw_data.loc[r, 'Correlation Data[1]':'Correlation Data[192]'])
            CDT.index=CD.index
            CDT_CD = pd.merge(CDT,CD,left_index=True,right_index=True).astype(float)
            CDT_CD.columns=np.arange(2)
            CDT_CD.index=np.arange(len(CDT_CD))
        
            # Fit
            popt, pcov = curve_fit(func_ICF,CDT_CD[0],CDT_CD[1], 
                                   p0=initial,bounds=(0, [1,1, np.inf, np.inf,1]))
            
            # Calculate R2 value
            yfitted=func_ICF(CDT,*popt)
            yfitted.index=np.arange(192)
            ss_res = np.sum((CDT_CD.iloc[0:192,1]-yfitted[1])**2)
            ss_tot = np.sum((CDT_CD.iloc[0:192,1]-np.mean(CDT_CD.iloc[0:192,1]))**2)
            R2 = 1 - (ss_res / ss_tot)
            
            # Store data
            popts[r]=popt
            sd[r]=np.sqrt(np.diag(pcov))
            r2.append(R2)
            
        except RuntimeError: # Process fitting error
            nans=np.zeros(5)
            nans[:]=np.nan
            popts[r]=nans
            sd[r]=nans
            print(' - Fitting failed',end="")

    # Final datat process
    r2_DF = pd.DataFrame(r2,columns=["R2"])
    popts.index=['sigma2', 'A', 'tauf', 'taus', 'beta']
    sd.index=['sd_sigma2', 'sd_A', 'sd_tauf', 'sd_taus', 'sd_beta']
    sd=sd.T
    popts=popts.T
    popts.append(r2)
    popts_sd = pd.merge(popts,sd,left_index=True,right_index=True)
    popts_sd_r2 = pd.merge(popts_sd,r2_DF,left_index=True,right_index=True)
        
    # Export parameters as csv
    exportname = filepath[:-5]+'_ICoFited_StretchedExponential.xlsx'
    popts_sd_r2.to_excel(exportname,index=False)
    print("\n  The result was exported as: "+exportname)


def PowerLaw_fit(filepath):  
    # Import the excel file
    raw_data = pd.read_excel(filepath)
    
    # names of each parameter corresponding to the paper
    def func_ICF(CDT, sigma2, A, tauf, taux, n):
        CD =  sigma2 * (A*np.e**(-CDT/tauf)+(1-A)*(1+(CDT/taux))**((n-1)/2))**2 * (tauf<taux)
        return CD

    # Initial value of fitting parameters
    initial=[0.85,0.4,100,300,0.6] #sigma2, A,tauf[us],taux[us],n
    print("  Initial values for the fitting")
    print("    Default:\n")
    print("\t sigma2 \t 0.85")
    print("\t A \t\t 0.4")
    print("\t tau_fast[µs] \t 100")
    print("\t tau_x[µs] \t 300")
    print("\t n \t\t 0.6\n")
    print("  Change from defaut? (y/n)")
    ini_change_determined = False
    while ini_change_determined == False:
        ini_change = input("    >> ")
        if ini_change == "n":
            ini_change_determined = True
        elif ini_change == "y":
            print("  Input new initial values:\n")
            initial[0] = float(input("\t sigma2 \t>> "))
            initial[1] = float(input("\t A \t\t>> "))
            initial[2] = float(input("\t tau_fast [µs] \t>> "))
            initial[3] = float(input("\t tau_x [µs] \t>> "))
            initial[4] = float(input("\t n \t\t>> "))
            ini_change_determined = True
        else:
            print("  Anser with \"y\" or \"n\"!")
    
    # Boundary condition
    print("\n  Boundary conditions:\n")
    print("\t 0 <\tsigma2\t< 1")
    print("\t 0 <\tA\t< 1")
    print("\t 0 <\ttau_fast")
    print("\t 0 <\ttau_x")
    print("\t 0 <\tn\t< 1")
    
    # Constrain
    print("\n  Constrain:\n")
    print("\t tau_fast < taus_x")
    
    # Fitted parameters and variances are collected in:
    popts=pd.DataFrame()
    sd=pd.DataFrame()
    r2=[]

    # delay time
    CDT = pd.DataFrame(raw_data.loc[1,'Correlation Delay Times[1] (µs)':'Correlation Delay Times[192] (µs)'])
    
    # Report the number of data to be fitted
    print("\n  ",len(raw_data),"data will be fitted.")
    
    # Fit
    for r in np.arange(0,len(raw_data)):
        print("\n  ",r+1,"/",len(raw_data),end="")
        try:
            # Make DataFrame containning CDT and CD
            CD = pd.DataFrame(raw_data.loc[r, 'Correlation Data[1]':'Correlation Data[192]'])
            CDT.index=CD.index
            CDT_CD = pd.merge(CDT,CD,left_index=True,right_index=True).astype(float)
            CDT_CD.columns=np.arange(2)
            CDT_CD.index=np.arange(len(CDT_CD))
        
            # Fit
            popt, pcov = curve_fit(func_ICF,CDT_CD[0],CDT_CD[1], 
                                   p0=initial,bounds=(0, [1,1, np.inf, np.inf,1]))
            
            # Calculate R2 value
            yfitted=func_ICF(CDT,*popt)
            yfitted.index=np.arange(192)
            ss_res = np.sum((CDT_CD.iloc[0:192,1]-yfitted[1])**2)
            ss_tot = np.sum((CDT_CD.iloc[0:192,1]-np.mean(CDT_CD.iloc[0:192,1]))**2)
            R2 = 1 - (ss_res / ss_tot)
            
            # Store data
            popts[r]=popt
            sd[r]=np.sqrt(np.diag(pcov))
            r2.append(R2)
            
        except RuntimeError: # Process fitting error
            nans=np.zeros(5)
            nans[:]=np.nan
            popts[r]=nans
            sd[r]=nans
            print(' - Fitting failed',end="")

    # Final datat process
    r2_DF = pd.DataFrame(r2,columns=["R2"])
    popts.index=['sigma2', 'A', 'tauf', 'taux', 'n']
    sd.index=['sd_sigma2', 'sd_A', 'sd_tauf', 'sd_taux', 'sd_n']
    sd=sd.T
    popts=popts.T
    popts.append(r2)
    popts_sd = pd.merge(popts,sd,left_index=True,right_index=True)
    popts_sd_r2 = pd.merge(popts_sd,r2_DF,left_index=True,right_index=True)
        
    # Export parameters as csv
    exportname = filepath[:-5]+'_ICoFited_Powerlaw.xlsx'
    popts_sd_r2.to_excel(exportname,index=False)
    print("\n  The result was exported as: "+exportname)





""" Main process """
### Ask which formula to use
print("""
Choose the formula and input the corresponding number (1 or 2).\n
  1 : Stretched exponential
       g2-1 = sigma^2 {Aexp[-tau/tau_f] + (1-A)exp[(tau/tau_s)^beta]}^2
  2 : Power-law
       g2-1 = sigma^2 {Aexp[-tau/tau_f] + (1-A)(1+(tau/tau_x))^((n-1)/2)}^2
""")
functype_determined = False
while functype_determined == False:
    functype = input("  >> ")
    if functype == str(1):
        print("  Formula: Stretched Exponential")
        functype_determined = True
    elif functype == str(2):
        print("  Formula: Power-law")
        functype_determined = True
    else:
        print("  Anser with \"1\" or \"2\"!")

### Choose excel file
files = glob.glob("*.xlsx")
print("""
Choose one file for fitting. File with the extension \".xlsx\" can only be detected.
""")
for i in range(len(files)):
    print("  "+str(i+1)+" : "+files[i]) # print the list of files
print()
file_determined = False    
while file_determined == False:
    filenum = int(input("  >> "))
    filenum = filenum-1 
    try:
        print("  File: "+files[int(filenum)])
        file_determined = True
    except:
        print("  Anser with integer in the list above!")

### Calculaions
print("""
Fitting Calculation
""")
if functype==str(1):
    StretchedExp_fit(filepath=files[int(filenum)])
elif functype==str(2):
    PowerLaw_fit(filepath=files[int(filenum)])
else:
    print("Unknown Error")

print("""
All process finished.
Copy and paste this window if you make a log of the calculation. 
""")

input("Press ENTER to close this window.")