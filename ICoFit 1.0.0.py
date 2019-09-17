""" ver. 1.0.0
    This is the first version of ICoFit. """
# Inport of necessary modules
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import wx

class FileDropTarget(wx.FileDropTarget):
    """ Drag & Drop Class """
    def __init__(self, window):
        wx.FileDropTarget.__init__(self)
        self.window = window

    def OnDropFiles(self, x, y, files):
        self.window.text_entry.SetLabel(files[0])

        return 0


class App(wx.Frame):
    """ GUI """
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(500, 400), style=wx.DEFAULT_FRAME_STYLE)

        # Panel
        p = wx.Panel(self, wx.ID_ANY)

        label = wx.StaticText(p, wx.ID_ANY, 'Drop the DLS data file here\n or input the path of the file below.\nExcel (.xlsx) is only allowed.', style=wx.SIMPLE_BORDER | wx.TE_CENTER)
        font = wx.Font(20, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        label.SetFont(font)
        label.SetBackgroundColour("#e0ffe0")

        # Set the drop target
        label.SetDropTarget(FileDropTarget(self))

        # Textbox
        self.text_entry = wx.TextCtrl(p, wx.ID_ANY)
        
        # Button
        btn = wx.Button(p, wx.ID_ANY, 'Fit')
        btn.Bind(wx.EVT_BUTTON, self.clicked)

        # Space to show fixed file path
        self.pathfixed = wx.StaticText(p, wx.ID_ANY,'')

        # Space to show process situation
        self.situ = wx.StaticText(p, wx.ID_ANY,'\n')

        # Layout
        layout = wx.BoxSizer(wx.VERTICAL)
        layout.Add(label, flag=wx.EXPAND | wx.ALL, border=10, proportion=1)
        layout.Add(self.text_entry, flag=wx.EXPAND | wx.ALL, border=10)
        layout.Add(btn, flag=wx.EXPAND |wx.ALL, border=10)
        layout.Add(self.pathfixed, flag=wx.EXPAND | wx.ALL, border=10)
        layout.Add(self.situ, flag=wx.EXPAND | wx.ALL, border=10)

        p.SetSizer(layout)

        self.Show()
    
    def clicked(self, event):
        # Extract the file path
        filepath = self.text_entry.GetValue()
        self.text_entry.Clear()
        self.pathfixed.SetLabel(filepath)
        self.situ.SetLabel('Calculating...\nThis might take minutes to be finished (Do not drop another file now.)')
        ICFfit1(filepath)
        self.situ.SetLabel('Omedetou! All calculations finished.\nYou can drop another file.')


class ICFfit1():
    ''' Fitting Class'''
    def __init__(self,filepath):
        
        # Inport excel
        raw_data = pd.read_excel(filepath)

        # names of each parameter corresponding to the paper
        def func_ICF(CDT, sigma2, A, tauf, taus, beta):
            CD =  sigma2 * (A*np.e**(-CDT/tauf)+(1-A)*np.e**(-(CDT/taus)**beta))**2
            return CD

        # Initial value of fitting parameters
        initial=[1,1,10**(-4),10**(-3),0.5] #sigma2, A,tauf,taus,beta
        popt=initial

        # Fitted parameters and variances are collected in:
        popts=pd.DataFrame()
        sd=pd.DataFrame()

        # time as sec. (not micro sec.)
        CDT = pd.DataFrame(raw_data.loc[1,'Correlation Delay Times[1] (µs)':'Correlation Delay Times[192] (µs)'])*10**(-6)

        # Fit
        for r in np.arange(0,len(raw_data)):
            try:
                # Make DataFrame containning CDT and CD
                CD = pd.DataFrame(raw_data.loc[r, 'Correlation Data[1]':'Correlation Data[192]'])
                CDT.index=CD.index
                CDT_CD = pd.merge(CDT,CD,left_index=True,right_index=True).astype(float)
                CDT_CD.columns=np.arange(2)
                CDT_CD.index=np.arange(len(CDT_CD))
                # Use the last fitted parameters as initial value for fitting
                initial=popt #sigma2, A,tauf,taus,beta
                # Fit
                popt, pcov = curve_fit(func_ICF,CDT_CD[0],CDT_CD[1],
                                       p0=initial,bounds=(0, [1,1, 10, 10,1]))
                popts[r]=popt
                sd[r]=np.sqrt(np.diag(pcov))
                # Report if fitting successed
                
            except RuntimeError: # Process fitting error
                nans=np.zeros(5)
                nans[:]=np.nan
                popts[r]=nans
                sd[r]=nans
                print('Fitting Error at row ',r,end='\t')

        popts.index=['sigma2', 'A', 'tauf', 'taus', 'beta']
        sd.index=['sd_sigma2', 'sd_A', 'sd_tauf', 'sd_taus', 'sd_beta']
        sd=sd.T
        popts=popts.T
        popts_sd = pd.merge(popts,sd,left_index=True,right_index=True)
        
        # Export parameters as csv
        exportname = filepath[:-5]+'_ICoFited.xlsx'
        popts_sd.to_excel(exportname,index=False)

app = wx.App()
App(None, -1, 'ICoFit 1.0.0 alpha')
app.MainLoop()
