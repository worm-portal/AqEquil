import os
import re
import sys
import shutil
import copy
import collections
import pickle

import warnings
from subprocess import Popen
import pkg_resources
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

def load(filename, messages=True):
    """
    Load a speciation file.

    Parameters
    ----------
    filename : str
        Name of the speciation file.

    Returns
    ----------
    An object of class `Speciation`.
    """

    with open(filename, 'rb') as handle:
        speciation = pickle.load(handle)
        if messages:
            print("Loaded '{}'".format(filename))
        return speciation


def convert_to_RVector(value, force_Rvec=True):
    
    """
    Convert a value or list into an R vector of the appropriate type.
    
    Parameters
    ----------
    value : numeric or str, or list of numeric or str
        Value to be converted.
    
    force_Rvec : bool, default True
        If `value` is not a list, force conversion into a R vector?
        False will return an int, float, or str if value is non-list.
        True will always return an R vector.
    
    Returns
    -------
    int, float, str, or an rpy2 R vector
        A value or R vector of an appropriate data type.
    """

    if not isinstance(value, list) and not force_Rvec:
        return value
    elif not isinstance(value, list) and force_Rvec:
        value = [value]
    else:
        pass

    if all(isinstance(x, bool) for x in value):
        return ro.BoolVector(value)
    elif all(isinstance(x, int) for x in value):
        return ro.IntVector(value)
    elif all(isinstance(x, float) or isinstance(x, int) for x in value):
        return ro.FloatVector(value)
    else:
        return ro.StrVector(value)


class Speciation(object):
    
    """
    Stores the output of a speciation calculation.
    
    Attributes
    ----------
    input : pd.Dataframe
        Pandas dataframe containing user-supplied sample chemistry data.
    
    mass_contribution : pd.Dataframe
        Pandas dataframe containing basis species contributions to mass balance
        of aqueous species.
    
    batch_3o : rpy2 ListVector
        An rpy2 ListVector (R object) containing speciation results, in case
        analysis in R is preferred.
    
    report : pd.Dataframe
        Pandas dataframe reporting major results of speciation calculation in
        across all samples.
    
    report_divs : rpy2 ListVector
        An rpy2 ListVector of column names within the different sections of the
        speciation report.
    
    sample_data : dict
        Dictionary with sample names as keys and speciation results as values.
    
    """
    
    def __init__(self, args):
        for k in args:
            setattr(self, k, args[k])

    def __getitem__(self, item):
         return getattr(self, item)
    
    @staticmethod
    def __unique(seq):
        """
        Provide a sequence, get a list of non-repeating elements in the same order.
        """
        seen = set()
        seen_add = seen.add
        return [x for x in seq if not (x in seen or seen_add(x))]

    
    def save(self, filename, messages=True):
        """
        Save the speciation as a '.speciation' file to your current working
        directory. This file can be loaded with `AqEquil.load(filename)`.
        
        Parameters
        ----------
        filename : str
            The desired name of the file.
            
        messages : str
            Print a message confirming the save?
        """
        
        if filename[-11:] != '.speciation':
            filename = filename + '.speciation'
        
        with open(filename, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)
            if messages:
                print("Saved as '{}'".format(filename))
    
    
    @staticmethod
    def __get_unit_info(subheader):
        
        unit_name_dict = {
            "pH" : ("", "pH"),
            "ppm" : ("", "ppm"),
            "ppb" : ("", "ppb"),
            "mg/L" : ("", "mg/L"),
            "degC" : ("temperature", "°C"),
            "log_molality" : ("log molality", "log(mol/kg)"),
            "Molality" : ("molality", "mol/kg"),
            "molality" : ("molality", "mol/kg"),
            "molal" : ("molality", "mol/kg"),
            "log_activity" : ("log activity", ""),
            "Log activity" : ("log activity", ""),
            "mg/kg.sol" : ("", "mg solute per kg solution"),
            "Alk., eq/kg.H2O" : ("alkalinity", "eq/kg"),
            "Alk., eq/L" : ("alkalinity", "eq/L"),
            "Alk., eq/kg.sol" : ("alkalinity", "eq/kg solution"),
            "Alk., mg/L CaCO3" : ("alkalinity", "mg/L CaCO3"),
            "Alk., mg/L HCO3-" : ("alkalinity", "mg/L HCO3-"),
            "pX" : ("-(log activity)", "-log(mol/kg)"),
            "activity" : ("activity", ""),
            "log_gamma" : ("log gamma", "log(kg/mol)"),
            "gamma" : ("gamma", "kg/mol"),
            "affinity_kcal" : ("affinity", "kcal/mol"),
            "%" : ("", "%"),
            "Eh_volts" : ("Eh", "volts"),
            "eq/kg.H2O" : ("charge", "eq/kg"),
            "logfO2" : ("", ""),
            "cal/kg.H2O" : ("energy supply", "cal/kg H2O"),
            "Log ion-H+ activity ratio" : ("Log ion-H+ activity ratio", ""),
        }
        
        out = unit_name_dict.get(subheader)
        
        return out[0], out[1]
    
    
    def lookup(self, col=None):
        
        """
        Look up desired columns in the speciation report.
        
        Parameters
        ----------
        col : str or list of str
            Leave blank to get a list of section names in the report:
            ```speciation.lookup()```
            Provide the name of a section to look up the names of columns in
            that section of the report:
            ```speciation.lookup("aq_distribution")```
            Provide a column name (or a list of column names) to retrieve the
            column from the report:
            ```speciation.lookup(["Temperature", "O2"])```
            
        Returns
        ----------
        Pandas dataframe or list of str
            If a column name (or list of column names) is provided, returns the
            speciation report with only the desired column(s). Otherwise returns
            a list of section names (if no arguments are provided), or a list of
            columns in a section (if a section name is provided).
        """
        
        if col==None and self.report_divs.named:
            return list(self.report_divs.names)
        
        if self.report_divs.named:
            if col in list(self.report_divs.names):
                return list(self.report_divs.rx2(col))
        
        if isinstance(col, str):
            col = [col]
        
        return self.report.iloc[:, self.report.columns.get_level_values(0).isin(set(col))]
    
    def __convert_aq_units_to_log_friendly(self, species, rows):

        col_data = self.lookup(species)
        
        col_data = col_data.loc[rows]
        
        if col_data.columns.get_level_values(1) == 'log_activity':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'activity'
        elif col_data.columns.get_level_values(1) == 'log_molality':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'molality'
        elif col_data.columns.get_level_values(1) == 'log_gamma':
            y = [10**float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = 'gamma'
        else:
            y = [float(s[0]) if s[0] != 'NA' else float("nan") for s in col_data.values.tolist()]
            out_unit = col_data.columns.get_level_values(1)[0]
        return y, out_unit
    
    
    def plot_mineral_saturation(self, sample_name, mineral_sat_type="affinity",
                                yrange=None,
                                colors=["blue", "orange"], bg_color="white",
                                save_as=None):
        """
        Vizualize mineral saturation states in a sample as a bar plot.
        
        Parameters
        ----------
        sample_name : str
            Name of the sample to plot.
        
        mineral_sat_type : str, default "affinity"
            Metric for mineral saturation state to plot. Can be "affinity" or
            "logQoverK".
            
        yrange : list of numeric, optional
            Sets the lower and upper limits of the y axis.
        
        colors : list of two str, default ["blue", "orange"]
            Sets the color of the bars representing supersaturated
            and undersaturated states, respectively.
        
        bg_color : str, default "white"
            Name of the Matplotlib color you wish to set as the panel
            background. A list of named colors can be found here:
            https://matplotlib.org/stable/gallery/color/named_colors.html
            
        save_as : str, optional
            Provide a filename to save this figure as a PNG.
        """
        
        if sample_name not in self.report.index:
            msg = ("Could not find '{}'".format(sample_name)+" among sample "
                   "names in the speciation report. Sample names include "
                   "{}".format(list(self.report.index)))
            raise Exception(msg)
        

        
        if isinstance(self.sample_data[sample_name].get('mineral_sat', None), pd.DataFrame):
            mineral_data = self.sample_data[sample_name]['mineral_sat'][mineral_sat_type].astype(float).sort_values(ascending=False)
            x = mineral_data.index
        else:
            msg = ("This sample does not have mineral saturation state data."
                   "To generate this data, ensure get_mineral_sat=True when "
                   "running speciate(), or ensure this sample has "
                   "mineral-forming basis species.")
            raise Exception(msg)
            
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        plt.xticks(rotation = 45, ha='right')
        
        pos_sat = [m if m >= 0 else float("nan") for m in mineral_data] # possibly: special list for m==0
        neg_sat = [m if m < 0 else float("nan") for m in mineral_data]
        
        barlist = [] # stores sets of bars in the bar chart so they can referenced for annotation
        for i, y_plot in enumerate([pos_sat, neg_sat]):
            
            if i == 0:
                color = colors[0]
            else:
                color = colors[1]
                
            bars = ax.bar(x, y_plot, tick_label=x, color=color)
            
            barlist.append(bars)
            
            if mineral_sat_type == "affinity":
                ylabel = 'affinity, kcal/mol'
            if mineral_sat_type == "logQoverK":
                ylabel = 'logQ/K'
        
        if yrange != None:
            plt.ylim(yrange[0], yrange[1])
        
        ax.set_facecolor(bg_color)
        plt.ylabel(ylabel)

        if save_as != None:
            if ".png" not in save_as[:-4]:
                save_as = save_as+".png"
            
            plt.savefig(save_as, dpi=300, bbox_inches="tight")
            print("Saved figure as {}".format(save_as))
        
        plt.show()
    

    def barplot(self, y, yrange=None, show_trace=True,
                show_legend=True, show_missing=True, legend_loc="best",
                colormap="viridis", bg_color="white", save_as=None):
        
        """
        Show a bar plot to vizualize one or more variables across all samples.
        
        Parameters
        ----------
        y : str or list of str
            Name (or list of names) of the variables to plot. Valid variables
            are columns in the speciation report.
       
        yrange : list of numeric, optional
            Sets the lower and upper limits of the y axis.
        
        show_trace : bool, default True
            Show asterisks for columns with numerical values but are too short
            to see clearly?
            
        show_legend : bool, default True
            Show a legend if there is more than one variable?
        
        show_missing : bool, default True
            Show samples that do not have bars?
        
        legend_loc : str or pair of float, default "best"
            Location of the legend on the plot. See
            https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html#matplotlib.axes.Axes.legend
        
        colormap : str, default "viridis"
            Name of the Matplotlib colormap to color the barplot. See
            https://matplotlib.org/stable/tutorials/colors/colormaps.html

        bg_color : str, default "white"
            Name of the Matplotlib color you wish to set as the panel
            background. A list of named colors can be found here:
            https://matplotlib.org/stable/gallery/color/named_colors.html
            
        save_as : str, optional
            Provide a filename to save this figure as a PNG.
        """
        
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        plt.xticks(rotation = 45, ha='right')

        if not isinstance(y, list):
            y = [y]
        
        y_cols = self.lookup(y)
        
        if not show_missing:
            y_cols = y_cols.dropna(how='all')
        
        x = y_cols.index # names of samples
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(y)-1)
        cmap = cm.__getattribute__(colormap)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        X = np.arange(len(x))
        
        barlist = [] # stores sets of bars in the bar chart so they can referenced for annotation
        for i, yi in enumerate(y):
            y_col = y_cols.iloc[:, y_cols.columns.get_level_values(0)==yi]
            
            try:
                subheader = y_col.columns.get_level_values(1)[0]
            except:
                msg = ("Could not find '{}' ".format(yi)+"in the speciation "
                       "report. Available variables include "
                      "{}".format(list(set(self.report.columns.get_level_values(0)))))
                raise Exception(msg)
            try:
                unit_type, unit = self.__get_unit_info(subheader)
            except:
                unit_type = ""
                unit = ""
            
            try:
                y_vals = [float(y0[0]) if y0[0] != 'NA' else float("nan") for y0 in y_col.values.tolist()]
            except:
                msg = ("One or more the values belonging to "
                       "'{}' are non-numeric and cannot be plotted.".format(y_col.columns.get_level_values(0)[0]))
                raise Exception(msg)
            
            if [abs(y0) for y0 in y_vals] != y_vals: # convert to bar-friendly units if possible
                if subheader in ["log_activity", "log_molality", "log_gamma"]:
                    y_plot, out_unit = self.__convert_aq_units_to_log_friendly(yi, rows=x)
                    unit_type, unit = self.__get_unit_info(out_unit)
                else:
                    y_plot = y_vals
            else:
                y_plot = y_vals
            
            if i == 0:
                subheader_previous = subheader
                unit_previous = unit
            if unit != unit_previous and i != 0:
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format(unit, yi_previous, unit_previous)+""
                       "Plotted variables must share units.")
                raise Exception(msg)
            elif "activity" in subheader.lower() and "molality" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("activity", yi_previous, "molality")+""
                       "Plotted variables must share units.")
                raise Exception(msg)
            elif "molality" in subheader.lower() and "activity" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("molality", yi_previous, "activity")+""
                       "Plotted variables must share units.")
                raise Exception(msg)
                
            yi_previous = copy.deepcopy(yi)
            unit_previous = copy.deepcopy(unit)
            subheader_previous = copy.deepcopy(subheader)
            
            if len(y) != 1:
                color = m.to_rgba(i)
            else:
                color = "black"
            
            bars = ax.bar(X+i*(1/(len(y)+1)), y_plot, tick_label=x, color=color, width=1/(len(y)+1))
            
            barlist.append(bars)

        max_bar_height = 0
        for bars in barlist:
            for p in bars.patches:
                max_bar_height = np.nanmax([max_bar_height, abs(p.get_height())])
                
        for i,bars in enumerate(barlist):
            for p in bars.patches:
                if show_trace and abs(p.get_height())/max_bar_height <= 0.009:
                    if len(y) != 1:
                        color = m.to_rgba(i)
                    else:
                        color = "black"
                
                    plt.annotate("*",
                                  (p.get_x() + p.get_width() / 2., p.get_height()),
                                  ha = 'center', va = 'center', xytext = (0, 10),
                                  color=color,
                                  weight='bold',
                                  fontsize=18,
                                  textcoords = 'offset points')
        
        if len(y) > 1:
            if unit != "":
                ylabel = "{} [{}]".format(unit_type, unit)
            else:
                ylabel = unit_type
            if show_legend:
                ax.legend(labels=y, loc=legend_loc)
        else:
            if 'pH' in y:
                ylabel = 'pH'
            elif 'Temperature' in y:
                ylabel = 'Temperature [°C]'
            else:
                if unit != "":
                    ylabel = "{} {} [{}]".format(y[0], unit_type, unit)
                else:
                    ylabel = "{} {}".format(y[0], unit_type)
        
        if yrange != None:
            plt.ylim(yrange[0], yrange[1])
        
        ax.set_facecolor(bg_color)
        plt.ylabel(ylabel)

        if save_as != None:
            if ".png" not in save_as[:-4]:
                save_as = save_as+".png"
            
            plt.savefig(save_as, dpi=300, bbox_inches="tight")
            print("Saved figure as {}".format(save_as))
        
        plt.show()
    
    
    def scatterplot(self, x="pH", y="Temperature", xrange=None, yrange=None,
                show_legend=True, legend_loc="best",
                colormap="viridis", bg_color="white", save_as=None):
        
        """
        Vizualize two or more sample variables with a scatterplot.
        
        Parameters
        ----------
        x, y : str, default for x is "pH", default for y is "Temperature"
            Names of the variables to plot against each other. Valid variables
            are columns in the speciation report. `y` can be a list of
            of variable names for a multi-series scatterplot.
       
        xrange, yrange : list of numeric, optional
            Sets the lower and upper limits of the x and y axis.
            
        show_legend : bool, default True
            Show a legend if there is more than one variable?
        
        legend_loc : str or pair of float, default "best"
            Location of the legend on the plot. See
            https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html#matplotlib.axes.Axes.legend
        
        colormap : str, default "viridis"
            Name of the Matplotlib colormap to color the scatterpoints. See
            https://matplotlib.org/stable/tutorials/colors/colormaps.html

        bg_color : str, default "white"
            Name of the Matplotlib color you wish to set as the panel
            background. A list of named colors can be found here:
            https://matplotlib.org/stable/gallery/color/named_colors.html
            
        save_as : str, optional
            Provide a filename to save this figure as a PNG.
        """

        if not isinstance(y, list):
            y = [y]
        
        if not isinstance(x, str):
            raise Exception("x must be a string.")
        
        x_col = self.lookup(x)
        
        try:
            xsubheader = x_col.columns.get_level_values(1)[0]
        except:
            msg = ("Could not find '{}' ".format(x)+"in the speciation "
                   "report. Available variables include "
                   "{}".format(list(set(self.report.columns.get_level_values(0)))))
            raise Exception(msg)
            
        try:
            x_plot = [float(x0[0]) if x0[0] != 'NA' else float("nan") for x0 in x_col.values.tolist()]
        except:
            msg = ("One or more the values belonging to "
                   "'{}' are non-numeric and cannot be plotted.".format(x_col.columns.get_level_values(0)[0]))
            raise Exception(msg)
        
        try:
            xunit_type, xunit = self.__get_unit_info(xsubheader)
        except:
            xunit_type = ""
            xunit = ""
        
        norm = matplotlib.colors.Normalize(vmin=0, vmax=len(y)-1)
        cmap = cm.__getattribute__(colormap)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)
        
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        
        for i, yi in enumerate(y):
            y_col = self.lookup(yi)
            
            try:
                subheader = y_col.columns.get_level_values(1)[0]
            except:
                msg = ("Could not find '{}' ".format(yi)+"in the speciation "
                       "report. Available variables include "
                      "{}".format(list(set(self.report.columns.get_level_values(0)))))
                raise Exception(msg)
            try:
                unit_type, unit = self.__get_unit_info(subheader)
            except:
                unit_type = ""
                unit = ""
            
            try:
                y_plot = [float(y0[0]) if y0[0] != 'NA' else float("nan") for y0 in y_col.values.tolist()]
            except:
                msg = ("One or more the values belonging to "
                       "'{}' are non-numeric and cannot be plotted.".format(y_col.columns.get_level_values(0)[0]))
                raise Exception(msg)
                
            if i == 0:
                subheader_previous = subheader
                unit_previous = unit
            if unit != unit_previous and i != 0:
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format(unit, yi_previous, unit_previous)+""
                       "Plotted variables must share units.")
                raise Exception(msg)
            elif "activity" in subheader.lower() and "molality" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("activity", yi_previous, "molality")+""
                       "Plotted variables must share units.")
                raise Exception(msg)
            elif "molality" in subheader.lower() and "activity" in subheader_previous.lower():
                msg = ("{} has a different unit of measurement ".format(yi)+""
                       "({}) than {} ({}). ".format("molality", yi_previous, "activity")+""
                       "Plotted variables must share units.")
                raise Exception(msg)
                
            yi_previous = copy.deepcopy(yi)
            unit_previous = copy.deepcopy(unit)
            subheader_previous = copy.deepcopy(subheader)
            
            if len(y) != 1:
                color = m.to_rgba(i)
            else:
                color = "black"
            
            plt.scatter(x_plot, y_plot, marker='o', color=color)

        if len(y) > 1:
            if unit != "":
                ylabel = "{} [{}]".format(unit_type, unit)
            else:
                ylabel = unit_type
            if show_legend:
                ax.legend(labels=y, loc=legend_loc)
        else:
            if 'pH' in y:
                ylabel = 'pH'
            elif 'Temperature' in y:
                ylabel = 'Temperature [°C]'
            else:
                if unit != "":
                    ylabel = "{} {} [{}]".format(y[0], unit_type, unit)
                else:
                    ylabel = "{} {}".format(y[0], unit_type)
        
        if x == 'pH':
            xlabel = 'pH'
        elif x == 'Temperature':
            xlabel = 'Temperature [°C]'
        else:
            if xunit != "":
                xlabel = "{} {} [{}]".format(x, xunit_type, xunit)
            else:
                xlabel = "{} {}".format(x, xunit_type)
        
        if xrange != None:
            plt.xlim(xrange[0], xrange[1])
        
        if yrange != None:
            plt.ylim(yrange[0], yrange[1])

        ax.set_facecolor(bg_color)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

        if save_as != None:
            if ".png" not in save_as[:-4]:
                save_as = save_as+".png"
            
            plt.savefig(save_as, dpi=300, bbox_inches="tight")
            print("Saved figure as {}".format(save_as))
        
        plt.show()
    
    
    def plot_mass_contribution(self, basis, sort_by=None, ascending=True,
                                     sort_y_by=None, width=0.9,
                                     legend_loc=(1.02, 0.5),
                                     save_as=None):
        
        """
        Plot basis species contributions to mass balance of aqueous species
        across all samples.
        
        Parameters
        ----------
        basis : str
            Name of the basis species.
            
        sort_by : str, optional
            Name of the variable used to sort samples. Variable names must be
            taken from the speciation report column names. No sorting is done by
            default.
        
        ascending : bool, default True
            Should sample sorting be in ascending order? Descending if False.
            Ignored unless `sort_by` is defined.
        
        sort_y_by : list of str, optional
            List of species names in the order that they should be stacked, from
            the bottom of the plot to the top.
        
        width : float, default 0.9
            Width of bars. No space between bars if width=1.0.
        
        legend_loc : str or pair of float, default (1.02, 0.5)
            Location of the legend on the plot. See
            https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.legend.html#matplotlib.axes.Axes.legend
        
        save_as : str, optional
            Provide a filename to save this figure as a PNG.
        """
        
        try:
            self.mass_contribution
        except:
            msg = ("Results for basis species contributions to aqueous mass "
                   "balance could not be found. Ensure that "
                   "get_mass_contribution = True when running speciate().")
            raise Exception(msg)
            
        if basis not in set(self.mass_contribution['basis']):
            msg = ("The basis species {} ".format(basis)+"could not be found "
                   "among available basis species: "
                   "{}".format(str(list(set(self.mass_contribution['basis'])))))
            raise Exception(msg)
            
        df_sp = copy.deepcopy(self.mass_contribution.loc[self.mass_contribution['basis'] == basis])
        
        if sort_by != None:
            if sort_by in self.report.columns.get_level_values(0):
                sort_col = self.lookup(sort_by)
                sort_by_unit = sort_col.columns.get_level_values(1)[0]
                sort_index = sort_col.sort_values([(sort_by, sort_by_unit)], ascending=ascending).index
                
                df_list = []
                for i in sort_index:
                    df_list.append(df_sp[df_sp['sample']==i])

                df_sp = pd.concat(df_list)
                
            else:
                msg = ("Could not find {}".format(sort_by)+" in the "
                       "speciation report. Available variables include "
                       "{}".format(list(self.report.columns.get_level_values(0))))
                raise Exception(msg)
        
        df_sp['percent'] = df_sp['percent'].astype(float)
        
        unique_species = self.__unique(df_sp["species"])
        
        if "Other" in unique_species:

            unique_species.append(unique_species.pop(unique_species.index("Other")))
        
        labels = self.__unique(df_sp["sample"])

        bottom = np.array([0]*len(labels))

        if sort_y_by != None:
            if isinstance(sort_y_by, list):
                if len(unique_species) == len(sort_y_by):
                    if len([s for s in unique_species if s in sort_y_by]) == len(unique_species) and len([s for s in sort_y_by if s in unique_species]) == len(unique_species):
                        unique_species = sort_y_by
                    else:
                        valid_needed = [s for s in unique_species if s not in sort_y_by]
                        invalid = [s for s in sort_y_by if s not in unique_species]
                        msg = ("sort_y_by is missing the following species: "
                               "{}".format(valid_needed)+" and was provided "
                               "these invalid species: {}".format(invalid))
                        raise Exception(msg)
                        
                elif len(sort_y_by) < len(unique_species):
                    msg = ("sort_y_by must have of all of the "
                           "following species: {}".format(unique_species)+". "
                           "You are missing {}".format([s for s in unique_species if s not in sort_y_by]))
                    raise Exception(msg)
                else:
                    msg = ("sort_y_by can only have the "
                           "following species: {}".format(unique_species)+".")
                    raise Exception(msg)
            else:
                raise Exception("sort_y_by must be a list of species names.")
                
        fig, ax = plt.subplots()
        
        for i,sp in enumerate(unique_species):
            percents = []
            for sample in labels:
                df_sample = df_sp[df_sp["sample"]==sample]
                try:
                    percent = df_sample[df_sample["species"]==sp]["percent"].iloc[0]
                    percents.append(percent)
                except:
                    percents.append(0.0)
            ax.bar(labels, percents, width, bottom=bottom, label=sp)
            bottom = bottom + np.array(percents)

        ax.set_ylabel('mole %')
        ax.set_title('Species accounting for mass balance of '+basis)
        plt.xticks(rotation = 45, ha='right')
        ax.legend(loc=legend_loc)
        
        if save_as != None:
            if ".png" not in save_as[:-4]:
                save_as = save_as+".png"
            
            plt.savefig(save_as, dpi=300, bbox_inches="tight")
            print("Saved figure as {}".format(save_as))
            
        plt.show()


class AqEquil():

    """
    Class containing functions to speciate aqueous water chemistry data using
    existing or custom thermodynamic datasets.
    
    Parameters
    ----------
    eq36da : str, defaults to path given by the environment variable EQ36DA
        Path to directory where data1 files are stored. 
        
    eq36co : str, defaults to path given by the environment variable EQ36CO
        Path to directory where EQ3 executables are stored.
    
    Attributes
    ----------
    eq36da : str
        Path to directory where data1 files are stored.
        
    eq36co : str
        Path to directory where EQ3 executables are stored.
        
    df_input_processed : pd.Dataframe
        Pandas dataframe containing user-supplied sample chemistry data that has
        been processed by `speciate`.
        
    out_dict : pd.Dataframe
        Pandas dataframe reporting results of last speciation calculation
        performed by `speciate`.
    
    verbose : int, 0, 1, or 2, default 1
        Level determining how many messages are returned during a
        calculation. 2 for all messages, 1 for errors or warnings only,
        0 for silent.
        
    """

    def __init__(self,
                 eq36da=os.environ.get('EQ36DA'),
                 eq36co=os.environ.get('EQ36CO')):

        self.eq36da = eq36da
        self.eq36co = eq36co
        self.df_input_processed = None
        self.out_dict = None
        self.verbose = 1

        os.environ['EQ36DA'] = self.eq36da  # set eq3 db directory
        os.environ['EQ36CO'] = self.eq36co  # set eq3 .exe directory
    
    
    @staticmethod
    def __file_exists(filename, ext='.csv'):
        """
        Check that a file exists and that it has the correct extension.
        Returns True if so, raises exception if not.
        """
        
        ext_dict = {
            ".csv" : "comma separated values (.csv)",
            ".txt" : "standard text (.txt)",
            ".rds" : "R Data (.rds)",
        }
        
        if ext in filename[-4:]:
            if os.path.exists(filename) and os.path.isfile(filename):
                return True
            else:
                err = "Cannot locate input file {}.".format(filename)
                raise Exception(err)
        else:
            err = ("Input file {}".format(filename) + " "
                "must be in {} format.".format(ext_dict[ext]))
            raise Exception(err)
        
        return False
    
    
    def _check_database_file(self, filename):
        
        """
        Check for problems in the thermodynamic database CSV.
        """
        
        # is the file a csv?
        self.__file_exists(filename)
        
        thermo_df = pd.read_csv(filename)
        
        # does this file have the proper headers?
        required_headers = ["name", "abbrv", "formula", "state",
                            "ref1", "ref2", "date", "E_units",
                            "G", "H", "S", "Cp", "V",
                            "a1.a", "a2.b", "a3.c", "a4.d", "c1.e", "c2.f",
                            "omega.lambda", "z.T",
                            "azero", "neutral_ion_type",
                            "dissrxn", "tag", "formula_ox"]
        
        missing_headers = []
        for header in required_headers:
            if header not in thermo_df.columns:
                missing_headers.append(header)
        if len(missing_headers) > 0:
            msg = ("The thermodynamic database file '{}'".format(filename)+" "
                   "is missing one or more required columns: "
                   "{}".format(", ".join(missing_headers))+". "
                   "Are these headers spelled correctly in the file?")
            raise Exception(msg)
        
        # does Cl-, O2(g), and O2 exist in the file?
        required_species = ["Cl-", "O2", "O2(g)"]
        missing_species = []
        for species in required_species:
            if species not in list(thermo_df["name"]):
                missing_species.append(species)
        if len(missing_species) > 0:
            msg = ("The thermodynamic database file '{}'".format(filename)+" "
                   "is missing required species:"
                   "{}".format(missing_species)+". Default thermodynamic values"
                   " will be used.")
            warnings.warn(msg)
        
        return

    
    def _check_sample_input_file(self, input_filename, exclude, db, custom_db,
                                       charge_balance_on, suppress_missing):
        """
        Check for problems in sample input file.
        """
        
        # does the input file exist? Is it a CSV?
        if self.__file_exists(input_filename):
            df_in = pd.read_csv(input_filename, header=None) # no headers for now so colname dupes can be checked
        else:
            raise Exception("_check_sample_input() error!")
        
        # are there any samples?
        if df_in.shape[0] <= 2:
            err_no_samples = ("The file {}".format(input_filename) + " "
                "must contain at least three rows: the "
                "first for column names, the second for column subheaders, "
                "followed by one or more rows for sample data.")
            raise Exception(err_no_samples)
        
        err_list = [] # for appending errors found in the sample input file
        
        # get header list
        col_list = list(df_in.iloc[0, 1:])
        
        # are there blank headers?
        if True in [isinstance(x, float) and x != x for x in col_list]:
            # isinstance(x, float) and x != x is a typesafe way to check for nan
            err_blank_header = ("One or more columns in the sample input "
                "file have blank headers. These might be empty columns. "
                "Only the first column may have a blank header. Remove any "
                "empty columns and/or give each header a name.")
            raise Exception(err_blank_header)
        
        # are there duplicate headers?
        dupe_cols = list(set([x for x in col_list if col_list.count(x) > 1]))
        if len(dupe_cols) > 0:
            err_dupe_cols = ("Duplicate column names are not allowed. "
                "Duplicate column names were found for:\n"
                "{}".format(str(dupe_cols)))
            err_list.append(err_dupe_cols)
        
        df_in.columns = df_in.iloc[0] # set column names
        df_in = df_in.drop(df_in.index[0], axis=0) # drop column name row
        df_in_headercheck = copy.deepcopy(df_in.iloc[:,1:]) # drop first column. Deepcopy slice because drop() doesn't work well with unnamed columns.
        
        # drop excluded headers
        for exc in exclude:
            if exc == df_in.columns[0]: # skip if 'sample' column is excluded
                continue
            try:
                df_in_headercheck = df_in_headercheck.drop(exc, axis=1) # drop excluded columns
            except:
                err_bad_exclude = ("Could not exclude the header '{}'".format(exc)+". "
                                   "This header could not be found in {}".format(input_filename)+"")
                err_list.append(err_bad_exclude)
        
        # get row list
        row_list = list(df_in.iloc[1:, 0])
        
        # are there blank rows?
        if True in [isinstance(x, float) and x != x for x in row_list]:
            # isinstance(x, float) and x != x is a typesafe way to check for nan
            err_blank_row = ("One or more rows in the sample input "
                "file have blank sample names. These might be empty rows. "
                "Remove any empty rows and/or give each sample a name. Sample "
                "names go in the first column.")
            raise Exception(err_blank_row)
            
        # are there duplicate rows?
        dupe_rows = list(set([x for x in row_list if row_list.count(x) > 1]))
        if len(dupe_rows) > 0:
            err_dupe_rows = ("Duplicate sample names are not allowed. "
                "Duplicate sample names were found for:\n"
                "{}".format(str(dupe_rows)))
            err_list.append(err_dupe_rows)
        
        # are column names valid entries in the database?
        if custom_db:
            data0_path = "data0." + db
        else:
            data0_path = self.eq36da + "/data0." + db
        if os.path.exists(data0_path) and os.path.isfile(data0_path):
            with open(data0_path) as data0:
                data0_lines = data0.readlines()
                start_index = [i+1 for i, s in enumerate(data0_lines) if '*  species name' in s]
                end_index = [i-1 for i, s in enumerate(data0_lines) if 'elements' in s]
                db_species = [i.split()[0] for i in data0_lines[start_index[0]:end_index[0]]]
                if charge_balance_on == 'pH':
                    err_charge_balance_on_pH = ("To balance charge on pH, use "
                        "charge_balance_on='H+'")
                    err_list.append(err_charge_balance_on_pH)
                elif charge_balance_on in ['Temperature', 'logfO2']:
                    err_charge_balance_invalid_type = ("Cannot balance charge "
                        "on {}.".format(charge_balance_on))
                    err_list.append(err_charge_balance_invalid_type)
                elif charge_balance_on != "none" and charge_balance_on not in list(set(df_in_headercheck.columns)):
                    err_charge_balance_invalid_sp = ("The species chosen for charge balance"
                        " '{}'".format(charge_balance_on)+""
                        " was not found among the headers of the sample input file.")
                    err_list.append(err_charge_balance_invalid_sp)
                for species in list(set(df_in_headercheck.columns)):
                    if species not in db_species and species not in ['Temperature', 'logfO2', 'pH']:
                        err_species_not_in_db = ("The species '{}'".format(species) + " "
                            "was not found in {}".format(data0_path) + ". "
                            "If the column contains data that should not be "
                            "included in the speciation calculation, add the "
                            "column name to the 'exclude' argument. Try "
                            "help(AqEquil.AqEquil.speciate) "
                            "for more information about 'exclude'.")
                        err_list.append(err_species_not_in_db)
                    elif species == 'pH':
                        err_species_pH = ("Please rename the 'pH' column in "
                            "the sample input file to 'H+' with the subheader "
                            "unit 'pH'.")
                        err_list.append(err_species_pH)
                    
        else:
            err_no_data0 = ("Could not locate {}.".format(data0_path) + " "
                "Unable to determine if column headers included in "
                "{} ".format(input_filename) + "match entries for species "
                "in the requested thermodynamic database '{}'.".format(db))
            err_list.append(err_no_data0)
        
        
        # are subheader units valid?
        subheaders = df_in_headercheck.iloc[0,]
        valid_subheaders = ["degC", "ppm", "ppb", "Suppressed", "Molality",
                            "Molarity", "mg/L", "mg/kg.sol", "Alk., eq/kg.H2O",
                            "Alk., eq/L", "Alk., eq/kg.sol", "Alk., mg/L CaCO3",
                            "Alk., mg/L HCO3-", "Log activity", "Log act combo",
                            "Log mean act", "pX", "pH", "pHCl", "pmH", "pmX",
                            "Hetero. equil.", "Homo. equil.", "Make non-basis",
                            "logfO2"]
        for i, subheader in enumerate(subheaders):
            if subheader not in valid_subheaders:
                err_valid_sub = ("The subheader '{}'".format(subheader) + " "
                    "for the column '{}'".format(df_in_headercheck.columns[i]) + " "
                    "is not recognized. Valid subheaders are {}".format(str(valid_subheaders)) + ". "
                    "If the column {}".format(df_in_headercheck.columns[i]) + " "
                    "contains data that is not meant for the "
                    "speciation calculation, add the column name "
                    "to the 'exclude' argument. Try help(AqEquil.AqEquil.speciate) "
                    "for more information about 'exclude'.")
                err_list.append(err_valid_sub)
            
        # is a 'Temperature' column present?
        if "Temperature" not in df_in_headercheck.columns and "Temperature" not in exclude:
            err_temp = ("The column 'Temperature' was not found in the input file. "
                "Please include a column with 'Temperature' in the first row, "
                "'degC' in the second row, and a temperature value for each "
                "sample in degrees Celsius.")
            err_list.append(err_temp)
        
        # raise exception that outlines all errors found
        if len(err_list) > 0:
            errs = "\n\n*".join(err_list)
            errs = ("The input file {}".format(input_filename)+" encountered"
                " errors:\n\n*" + errs)
            raise Exception(errs)
        
        return
        
        
    def __clear_eqpt_extra_output(self):
        
        """
        Deletes all EQPT output except data1.
        """
        
        if os.path.exists("eqpt_log.txt") and os.path.isfile("eqpt_log.txt"):
            os.remove("eqpt_log.txt")
        if os.path.exists("data1f.txt") and os.path.isfile("data1f.txt"):
            os.remove("data1f.txt")
        if os.path.exists("slist.txt") and os.path.isfile("slist.txt"):
            os.remove("slist.txt")

            
    def runeqpt(self, db, extra_eqpt_output=False):
        
        """
        Convert a data0 into a data1 file with EQPT.
        
        Parameters
        ----------
        db : str
            Three letter code of database.
        
        extra_eqpt_output : bool, default False
            Keep additional output files from EQPT? These files include
            eqpt_log.txt, data1f.txt, and slist.txt.
        """

        if os.path.exists("data0."+db) and os.path.isfile("data0."+db):
            pass
        else:
            raise Exception("Error: could not locate custom database",
                            "data0.{} in {}.".format(db, os.getcwd()))

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            os.remove("data1."+db)

        self.__clear_eqpt_extra_output()

        os.environ['EQ36DA'] = os.getcwd()

        args = ['/bin/csh', self.eq36co+'/runeqpt', db]

        try:
            self.__run_script_and_wait(args) # run EQPT
        except:
            os.environ['EQ36DA'] = self.eq36da
            raise Exception(
                "Error: EQPT failed to run on {}.".format("data0."+db))

        if os.path.exists("data1") and os.path.isfile("data1"):
            os.rename("data1", "data1."+db)
        if os.path.exists("output") and os.path.isfile("output"):
            os.rename("output", "eqpt_log.txt")
        if os.path.exists("data1f") and os.path.isfile("data1f"):
            os.rename("data1f", "data1f.txt")
        if os.path.exists("slist") and os.path.isfile("slist"):
            os.rename("slist", "slist.txt")

        if os.path.exists("data1."+db) and os.path.isfile("data1."+db):
            if self.verbose > 0:
                print("Successfully created a data1."+db+" from data0."+db)
        else:
            msg = ("EQPT could not create data1."+db+" from "
                   "data0."+db+". Check eqpt_log.txt for details.")
            raise Exception(msg)

        if not extra_eqpt_output:
            self.__clear_eqpt_extra_output()

        os.environ['EQ36DA'] = self.eq36da  # reset default EQ36 db path

        
    def runeq3(self, filename_3i, db,
               samplename=None,
               path_3i=os.getcwd(),
               path_3o=os.getcwd(),
               path_3p=os.getcwd()):
        
        """
        Call EQ3 on a .3i input file.
        
        Parameters
        ----------
        filename_3i : str
            Name of 3i input file.
        
        db : str
            Three letter code of database.
        
        path_3i : path str, default current working directory
            Path of .3i input files.
            
        path_3o : path str, default current working directory
            Path of .3o output files.
        
        path_3p : path str, default current working directory
            Path of .3p pickup files.
        """

        # get current working dir
        cwd = os.getcwd()
        
        if samplename == None:
            samplename = filename_3i[:-3]
        
        if self.verbose > 0:
            print('Using ' + db + ' to speciate ' + samplename)
        os.chdir(path_3i)  # step into 3i folder
        args = ['/bin/csh', self.eq36co+'/runeq3', db, filename_3i]

        self.__run_script_and_wait(args) # run EQ3

        # restore working dir
        os.chdir(cwd)

        filename_3o = filename_3i[:-1] + 'o'
        filename_3p = filename_3i[:-1] + 'p'

        try:
            # rename output
            os.rename(path_3i + '/output', path_3i + "/" + filename_3o)
        except:
            if self.verbose > 0:
                print('Error: EQ3 failed to produce output for ' + filename_3i)

        try:
            # move output
            shutil.move(path_3i + "/" + filename_3o,
                        path_3o + "/" + filename_3o)
        except:
            if self.verbose > 0:
                print('Error: Could not move', filename_3o, "to", path_3o)

        try:
            # rename pickup
            os.rename(path_3i + '/pickup', path_3i + "/" + filename_3p)
            move_pickup = True
        except:
            if self.verbose > 0:
                print('Error: EQ3 failed to produce a pickup file for ' + filename_3i)
            move_pickup = False
        
        if move_pickup:
            try:
                # move pickup
                shutil.move(path_3i + "/" + filename_3p,
                            path_3p + "/" + filename_3p)
            except:
                if self.verbose > 0:
                    print('Error: Could not move', filename_3p, "to", path_3p)

                    
    def runeq6(self, filename_6i, db,
               samplename=None,
               path_6i=os.getcwd(),
               path_6o=os.getcwd(),
               path_6p=os.getcwd()):
        
        """
        Call EQ6 on a .6i input file.
        
        Parameters
        ----------
        filename_6i : str
            Name of 6i input file.
        
        db : str
            Three letter code of database.
        
        path_6i : path str, default current working directory
            Path of .6i input files.
            
        path_6o : path str, default current working directory
            Path of .6o output files.
        
        path_6p : path str, default current working directory
            Path of .6p pickup files.
        """

        # get current working dir
        cwd = os.getcwd()
        
        if samplename == None:
            samplename = filename_6i[:-3]
        
        if self.verbose > 0:
            print('Using ' + db + ' to speciate ' + samplename)
        os.chdir(path_6i)  # step into 6i folder
        args = ['/bin/csh', self.eq36co+'/runeq6', db, filename_6i]

        self.__run_script_and_wait(args) # run EQ6

        # restore working dir
        os.chdir(cwd)

        filename_6o = filename_6i[:-1] + 'o'
        filename_6p = filename_6i[:-1] + 'p'

        try:
            # rename output
            os.rename(path_6i + '/output', path_6i + "/" + filename_6o)
        except:
            if self.verbose > 0:
                print('Error: EQ6 failed to produce output for ' + filename_6i)

        try:
            # move output
            shutil.move(path_6i + "/" + filename_6o,
                        path_6o + "/" + filename_6o)
        except:
            if self.verbose > 0:
                print('Error: Could not move', filename_6o, "to", path_6o)

        try:
            # rename pickup
            os.rename(path_6i + '/pickup', path_6i + "/" + filename_6p)
            move_pickup = True
        except:
            if self.verbose > 0:
                print('Error: EQ6 failed to produce a pickup file for ' + filename_6i)
            move_pickup = False
        
        if move_pickup:
            try:
                # move pickup
                shutil.move(path_6i + "/" + filename_6p,
                            path_6p + "/" + filename_6p)
            except:
                if self.verbose > 0:
                    print('Error: Could not move', filename_6p, "to", path_6p)
                    
                    
    def __mk_check_del_directory(self, path):
        
        """
        Checks for the dir being created. If it is already present, delete it
        before recreating it.
        """
        
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            shutil.rmtree(path)
            os.makedirs(path)

            
    def __read_inputs(self, file_type, location):
        
        """
        Finds all files of a filetype in all downstream folders.
        """
        
        file_name = []  # file names
        file_list = []  # file names with paths
        for root, dirs, files in os.walk(location):
            for file in files:
                if file.endswith(file_type):
                    if "-checkpoint" not in file:
                        file_name.append(file)
                        file_list.append(os.path.join(root, file))
        return file_name, file_list

    
    def __run_script_and_wait(self, args):
        
        """
        Runs shell commands.
        """
        
        with open(os.devnull, 'w') as fp:  # devnull supresses written output
            Popen(args, stdout=fp).wait()

            
    def _delete_rxn_folders(self):
        
        """
        Deletes folders storing raw EQ3/6 input and output.
        """
        
        if os.path.exists('rxn_3i') and os.path.isdir('rxn_3i'):
            shutil.rmtree('rxn_3i')
        if os.path.exists('rxn_3o') and os.path.isdir('rxn_3o'):
            shutil.rmtree('rxn_3o')
        if os.path.exists('rxn_3p') and os.path.isdir('rxn_3p'):
            shutil.rmtree('rxn_3p')
        if os.path.exists('rxn_6i') and os.path.isdir('rxn_6i'):
            shutil.rmtree('rxn_6i')
        if os.path.exists('rxn_6o') and os.path.isdir('rxn_6o'):
            shutil.rmtree('rxn_6o')
        if os.path.exists('rxn_6p') and os.path.isdir('rxn_6p'):
            shutil.rmtree('rxn_6p')
            

    def speciate(self,
                 input_filename,
                 db="wrm",
                 redox_flag=0,
                 redox_aux="Fe+3",
                 default_logfO2=-6,
                 exclude=[],
                 suppress=[],
                 charge_balance_on="none",
                 suppress_missing=True,
                 verbose=1,
                 report_filename=None,
                 get_aq_dist=True,
                 aq_dist_type="log_activity",
                 get_mass_contribution=True,
                 mass_contribution_other=True,
                 get_mineral_sat=True,
                 mineral_sat_type="affinity",
                 get_redox=True,
                 redox_type="Eh",
                 get_ion_activity_ratios=True,
                 get_affinity_energy=False,
                 rxn_filename=None,
                 not_limiting=["H+", "OH-", "H2O"],
                 get_charge_balance=True,
                 custom_db=False,
                 extra_eqpt_output=False,
                 batch_3o_filename=None,
                 delete_generated_folders=False):
        
        """
        Calculate the equilibrium distribution of chemical species in solution.
        Additionally, calculate chemical affinities and energy supplies for
        user-specified reactions.
        
        Parameters
        ----------
        input_filename : str
            User-supplied utf8-encoded comma separated value (csv) file
            containing sample data intended for speciation. The file must
            follow this format:
            
            - the first row is a header row that must contain the names of the
              species to be included in the speciation calculation. There
              cannot be duplicate headers.
            - the second row must contain subheaders for each species in the
              header row. These subheaders must be taken from the following:
              
                    degC
                    ppm
                    ppb
                    Suppressed
                    Molality
                    Molarity
                    mg/L
                    mg/kg.sol
                    Alk., eq/kg.H2O
                    Alk., eq/L
                    Alk., eq/kg.sol
                    Alk., mg/L CaCO3
                    Alk., mg/L HCO3-
                    Log activity
                    Log act combo
                    Log mean act
                    pX
                    pH
                    pHCl
                    pmH
                    pmX
                    Hetero. equil.
                    Homo. equil.
                    Make non-basis
                    
            - 'Temperature' must be included as a header, with 'degC' as its
              subheader.
            - The first column must contain sample names. There cannot be
              duplicate sample names.
        
        db : three letter str, default "wrm"
            Three letter file extension for the desired thermodynamic database.
            If `custom_db` is False, this database must be named data1.xyz
            (where xyz is your desired three letter extension) and located
            in the EQ3/6 'EQ36DA' path. Otherwise, the database must be named
            data0.xyz and located in your current working directory. Note that
            data1 files are already compiled by EQPT, while data0 files will be
            automatically compiled for you if `custom_db` is True.
        
        redox_flag : int, default 0
            Values corresponding to redox options in the EQ3/6 v8.0 software.
            For more information see the 'Redox Option' section of the EQ3/6
            version 8.0 software user's manual. Set sample redox state with the
            following options:
            
            * -3 for O2(g)
            * -2 for pe (in pe units)
            * -1 for Eh (volts)
            *  0 for logfO2 (log bars), or dissolved O2 (see below)
            *  1 for defining a redox couple (see `redox_aux`)
             
            Note that if you are importing water chemistry data from a
            spreadsheet, a column must be supplied with data that corresponds to
            the redox option you chose. The column name is important:
            
            * -3 must have a column named: O2(g)
            * -2 must have a column named: pe
            * -1 must have a column named: Eh
            *  0 must have a column named: logfO2
            *  1 must have a column corresponding to the auxilliary basis species
              selected to form a redox couple with its linked species (see
              `redox_aux`). For example, the redox couple Fe+2/Fe+3 would need
              a column named: Fe+3
            
            If an appropriate header or redox data cannot be found to define
            redox state, `default_logfO2` is used to set sample logfO2.
            
            There is a special case where dissolved oxygen can be used to impose
            sample redox state if `redox_flag` is set to 0 and a column named
            logfO2 does not appear in the sample data sheet. If there is a
            column corresponding to dissolved oxygen measurements, logfO2 is
            calculated from the equilibrium reaction O2(aq) = O2(g) at the
            temperature and pressure of the sample using the revised Helgeson-
            Kirkham-Flowers (HKF) equation of state (JC Tanger IV and HC
            Helgeson, Am. J. Sci., 1988, 288, 19).
        
        redox_aux : default "Fe+3", optional
            Ignored unless `redox_flag` equals 1. Name of the auxilliary species
            whose reaction links it to a basis species (or another auxilliary
            species) such that they form a redox couple that controls sample
            fO2. For instance, Fe+3 is linked to Fe+2 in many supporting data
            files, so selecting `redox_flag` = 1 and `redox_aux` = "Fe+3" will
            set sample fO2 based on the Fe+2/Fe+3 redox couple.
        
        default_logfO2 : float, default -6
            Default value for sample logfO2 in case redox data cannot be found
            in the user-supplied sample spreadsheet.
        
        exclude : list of str, default []
            Names of columns in the user-supplied sample spreadsheet that should
            not be considered aqueous species. Useful for excluding columns
            containing sample metatadata, such as "Year" and "Location".
            
        suppress : list of str, default []
            Names of chemical species that will be prevented from forming in the
            speciation calculation.
            
        charge_balance_on : str, default "none"
            If "none", will not balance electrical charge between cations and
            anions in the speciation calculation. If a name of a species is
            supplied instead, the activity of that species will be allowed to
            change until charge balance is obtained. For example,
            charge_balance_on = "H+" will calculate what pH a sample must have
            to have zero net charge.
        
        suppress_missing : bool, default True
            Suppress the formation of an aqueous species if it is missing a
            value in the user-supplied sample spreadsheet?
        
        verbose : int, 0, 1, or 2, default 1
            Level determining how many messages are returned during a
            calculation. 2 for all messages, 1 for errors or warnings only,
            0 for silent.
            
        report_filename : str, optional
            Name of the comma separated values (csv) report file generated when
            the calculation is complete. If this argument is not defined, a
            report file is not generated.
            
        get_aq_dist : bool, default True
            Calculate distributions of aqueous species?
        
        aq_dist_type : str, default "log_activity"
            Desired units of measurement for reported distributions of aqueous
            species. Can be "molality", "log_molality", "log_gamma", or
            "log_activity". Ignored if `get_aq_dist` is False.
        
        get_mass_contribution : bool, default True
            Calculate basis species contributions to mass balance of aqueous
            species?
        
        mass_contribution_other : bool, default True
            Include an "other" species for the sake of summing percents of basis
            species contributions to 100%? Ignored if `get_mass_contribution` is
            False.
        
        get_mineral_sat : bool, default True
            Calculate saturation states of pure solids?
        
        mineral_sat_type : str, default "affinity"
            Desired units of measurement for reported saturation states of pure
            solids. Can be "logQoverK" or "affinity". Ignored if
            `get_mineral_sat` is False.
        
        get_redox : bool, default True
            Calculate potentials of redox couples?
            
        redox_type : str, default "Eh"
            Desired units of measurement for reported redox potentials. Can be
            "Eh", "pe", "logfO2", or "Ah". Ignored if `get_redox` is False.
        
        get_ion_activity_ratios : bool, default True
            Calculate ion/H+ activity ratios and neutral species activities?
        
        get_affinity_energy : bool, default False
            Calculate affinities and energy supplies of reactions listed in a
            separate user-supplied file?
        
        rxn_filename : str, optional
            Name of .txt file containing reactions used to calculate affinities
            and energy supplies. Ignored if `get_affinity_energy` is False.
        
        not_limiting : list, default ["H+", "OH-", "H2O"]
            List containing names of species that are not considered limiting
            when calculating energy supplies. Ignored if `get_affinity_energy`
            is False.
        
        get_charge_balance : bool, default True
            Calculate charge balance and ionic strength?
            
        custom_db : bool, default False
            Is the database defined by `db` a custom user-supplied database? If
            this is set to True, searches for a data0.xyz file in the current
            working directory, where 'xyz' corresponds to the three letter code
            assigned to `db`. This data0 file is automatically converted into a
            machine-readable file called data1 by software called EQPT. This
            data1 file is then used in speciation calculations.
        
        extra_eqpt_output : bool, default False
            Keep additional output files created by EQPT (see `custom_db`)?
            Ignored if `custom_db` is False.
        
        batch_3o_filename : str, optional
            Name of rds (R object) file exported after the speciation
            calculation? No file will be generated if this argument is not
            defined.
            
        delete_generated_folders : bool, default False
            Delete the 'rxn_3i', 'rxn_3o', and 'rxn_3p' folders containing raw
            EQ3NR input, output, and pickup files once the speciation
            calculation is complete?
        
        Returns
        -------
        speciation : object of class Speciation
            Contains the results of the speciation calculation.
        
        """
        
        self.verbose = verbose
        
        # check input sample file for errors
        self._check_sample_input_file(input_filename, exclude, db, custom_db,
                                      charge_balance_on, suppress_missing)
        
        # handle batch_3o naming
        if batch_3o_filename != None:
            if ".rds" in batch_3o_filename[-4:]:
                batch_3o_filename = batch_3o_filename
            else:
                batch_3o_filename = "batch_3o_{}.rds".format(db)
        else:
            batch_3o_filename = ro.r("NULL")

        if custom_db:
            # EQ3/6 cannot handle spaces in the 'EQ36DA' path name.
            if " " in os.getcwd():
                msg = ("Error: the path to the custom database "
                    "cannot contain spaces. The current path "
                    "is: [ " + os.getcwd() + " ]. Remove or "
                    "replace spaces in folder names for this "
                    "feature. Example: [ " + os.getcwd().replace(" ", "-") + " ].")
                raise Exception(msg)

            self.runeqpt(db, extra_eqpt_output)
            os.environ['EQ36DA'] = os.getcwd()

        if get_affinity_energy:
            if rxn_filename == None:
                err = ("A get_affinity_energy was set to True but a reaction "
                       "file was not specified.")
                raise Exception(err)
            elif self.__file_exists(rxn_filename, '.txt'):
                pass
        else:
            rxn_filename = ""

        # preprocess for EQ3 using R scripts
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_prescript = pkg_resources.resource_string(
                __name__, 'preprocess_for_EQ3.r').decode("utf-8")
            ro.r(r_prescript)
            df_input_processed = ro.r.preprocess(input_filename=input_filename,
                                                 exclude=convert_to_RVector(
                                                     exclude),
                                                 redox_flag=redox_flag,
                                                 default_logfO2=default_logfO2,
                                                 charge_balance_on=charge_balance_on,
                                                 suppress_missing=suppress_missing,
                                                 suppress=convert_to_RVector(
                                                     suppress),
                                                 verbose=self.verbose)

        for warning in w:
            print(warning.message)

        self.df_input_processed = pandas2ri.ri2py_dataframe(df_input_processed)

        # run EQ3 on each input file
        cwd = os.getcwd()

        self.__mk_check_del_directory('rxn_3o')
        self.__mk_check_del_directory('rxn_3p')
        files_3i, files_3i_paths = self.__read_inputs('3i', 'rxn_3i')

        input_dir = cwd + "/rxn_3i/"
        output_dir = cwd + "/rxn_3o/"
        pickup_dir = cwd + "/rxn_3p/"
        
        for file in files_3i:
            samplename = self.df_input_processed.loc[file[:-3], "Sample"]
            self.runeq3(filename_3i=file, db=db, samplename=samplename,
                        path_3i=input_dir, path_3o=output_dir,
                        path_3p=pickup_dir)

        if custom_db:
            os.environ['EQ36DA'] = self.eq36da

        files_3o = [file+".3o" for file in self.df_input_processed.index]
        
        df_input_processed_names = convert_to_RVector(list(self.df_input_processed.columns))
        
        # mine output
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_3o_mine = pkg_resources.resource_string(
                __name__, '3o_mine.r').decode("utf-8")
            ro.r(r_3o_mine)
            batch_3o = ro.r.main_3o_mine(
                files_3o=convert_to_RVector(files_3o),
                input_filename=input_filename,
                rxn_filename=rxn_filename,
                get_aq_dist=get_aq_dist,
                aq_dist_type=aq_dist_type,
                get_mass_contribution=get_mass_contribution,
                mass_contribution_other=mass_contribution_other,
                get_mineral_sat=get_mineral_sat,
                mineral_sat_type=mineral_sat_type,
                get_redox=get_redox,
                redox_type=redox_type,
                get_charge_balance=get_charge_balance,
                get_ion_activity_ratios=get_ion_activity_ratios,
                get_affinity_energy=get_affinity_energy,
                not_limiting=convert_to_RVector(not_limiting),
                batch_3o_filename=batch_3o_filename,
                df_input_processed=pandas2ri.py2ri(self.df_input_processed),
                # Needed for keeping symbols in column names after porting
                #   df_input_processed in the line above. Some kind of check.names
                #   option for pandas2ri.py2ri would be nice. Workaround:
                df_input_processed_names=df_input_processed_names,
                verbose=self.verbose,
            )
        for warning in w:
            print(warning.message)
        
        if get_mass_contribution:
            mass_contribution = pandas2ri.ri2py_dataframe(batch_3o.rx2('mass_contribution'))
        df_report = pandas2ri.ri2py_dataframe(batch_3o.rx2('report'))
        df_input = pandas2ri.ri2py_dataframe(batch_3o.rx2('input'))
        report_divs = batch_3o.rx2('report_divs')

        input_cols = list(report_divs.rx2('input'))
        df_input = df_report.loc[:, input_cols]

        # handle headers and subheaders of input section
        headers = [col.split("_")[0] for col in list(df_input.columns)]
        headers = ["pH" if header == "H+" else header for header in headers]
        headers = [header+"_(input)" if header not in ["Temperature", "pH", "logfO2"]+exclude else header for header in headers]
        subheaders = [subheader[1] if len(subheader) > 1 else "" for subheader in [
            col.split("_") for col in list(df_input.columns)]]
        multicolumns = pd.MultiIndex.from_arrays(
            [headers, subheaders], names=['Sample', ''])
        df_input.columns = multicolumns

        df_join = df_input

        if get_aq_dist:
            aq_distribution_cols = list(report_divs.rx2('aq_distribution'))
            df_aq_distribution = df_report.loc[:, aq_distribution_cols]
            df_aq_distribution = df_aq_distribution.apply(pd.to_numeric, errors='coerce')

            # handle headers of aq_distribution section
            headers = df_aq_distribution.columns
            subheaders = [aq_dist_type]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_aq_distribution.columns = multicolumns
            df_join = df_join.join(df_aq_distribution)

        if get_mineral_sat:
            mineral_sat_cols = list(report_divs.rx2('mineral_sat'))
            df_mineral_sat = df_report.loc[:, mineral_sat_cols]
            df_mineral_sat = df_mineral_sat.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_mineral_sat section
            if mineral_sat_type == "affinity":
                mineral_sat_unit = "affinity_kcal"
            elif mineral_sat_type == "logQoverK":
                mineral_sat_unit = "logQ/K"
            else:
                raise Exception(
                    "mineral_sat_type must be either 'affinity' or 'logQoverK'")

            headers = df_mineral_sat.columns
            subheaders = [mineral_sat_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_mineral_sat.columns = multicolumns
            df_join = df_join.join(df_mineral_sat)

        if get_redox:
            redox_cols = list(report_divs.rx2('redox'))
            df_redox = df_report.loc[:, redox_cols]
            df_redox = df_redox.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_redox section
            if redox_type == "Eh":
                redox_unit = "Eh_volts"
            elif redox_type == "pe":
                redox_unit = "pe"
            elif redox_type == "logfO2":
                redox_unit = "logfO2"
            elif redox_type == "Ah":
                redox_unit = "Ah_kcal"
            else:
                raise Exception(
                    "redox_type must be either 'Eh', 'pe', 'logfO2', or 'Ah'")

            headers = df_redox.columns
            subheaders = [redox_unit]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_redox.columns = multicolumns
            df_join = df_join.join(df_redox)

        if get_charge_balance:
            charge_balance_cols = list(report_divs.rx2('charge_balance'))
            df_charge_balance = df_report.loc[:, charge_balance_cols]
            df_charge_balance = df_charge_balance.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_charge_balance section
            headers = df_charge_balance.columns
            subheaders = ["%"]*2 + ['eq/kg.H2O', 'molality'] + \
                ['eq/kg.H2O']*4 + ['molality']
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_charge_balance.columns = multicolumns
            df_join = df_join.join(df_charge_balance)
            
        if get_ion_activity_ratios:
            ion_activity_ratio_cols = list(report_divs.rx2('ion_activity_ratios'))
            df_ion_activity_ratios = df_report.loc[:, ion_activity_ratio_cols]
            df_ion_activity_ratios = df_ion_activity_ratios.apply(pd.to_numeric, errors='coerce')
            
            # handle headers of df_ion_activity_ratios section
            headers = df_ion_activity_ratios.columns
            subheaders = ["Log ion-H+ activity ratio"]*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_ion_activity_ratios.columns = multicolumns
            df_join = df_join.join(df_ion_activity_ratios)

        if get_affinity_energy:
            affinity_cols = list(report_divs.rx2('affinity'))
            energy_cols = list(report_divs.rx2('energy'))
            df_affinity = df_report.loc[:, affinity_cols]
            df_energy = df_report.loc[:, energy_cols]
            df_affinity = df_affinity.apply(pd.to_numeric, errors='coerce')
            df_energy = df_energy.apply(pd.to_numeric, errors='coerce')

            # handle headers of df_affinity section
            headers = df_affinity.columns
            subheaders = ['cal/mol e-']*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_affinity.columns = multicolumns

            # handle headers of df_energy section
            headers = df_energy.columns
            subheaders = ['cal/kg.H2O']*len(headers)
            multicolumns = pd.MultiIndex.from_arrays(
                [headers, subheaders], names=['Sample', ''])
            df_energy.columns = multicolumns
            df_join = df_join.join(df_affinity)
            df_join = df_join.join(df_energy)

        out_dict = {'sample_data': {},
                    'report': df_join,
                    'input': df_input, 'report_divs': report_divs}
        
        if get_mass_contribution:
            out_dict['mass_contribution'] = mass_contribution

        sample_data = batch_3o.rx2('sample_data')

        # assemble sample data
        for i, sample in enumerate(sample_data):
            dict_sample_data = {
                "filename": str(sample.rx2('filename')[0]),
                "name": str(sample.rx2('name')[0]),
                "temperature": float(sample.rx2('temperature')[0]),
                "pressure": float(sample.rx2('pressure')[0]),
                "logact_H2O": float(sample.rx2('logact_H2O')[0]),
                "H2O_density": float(sample.rx2('H2O_density')[0]),
                "H2O_molality": float(sample.rx2('H2O_molality')[0]),
                "H2O_log_molality": float(sample.rx2('H2O_log_molality')[0]),
                }

            if get_aq_dist:
                sample_aq_dist = pandas2ri.ri2py_dataframe(sample.rx2('aq_distribution'))
                sample_aq_dist = sample_aq_dist.apply(pd.to_numeric, errors='coerce')
                dict_sample_data.update({"aq_distribution": sample_aq_dist})

            if get_mass_contribution:
                sample_mass_contribution = mass_contribution[mass_contribution["sample"] == sample.rx2('name')[0]]
                dict_sample_data.update(
                    {"mass_contribution": sample_mass_contribution})

            if get_mineral_sat:
                dict_sample_data.update(
                    {"mineral_sat": pandas2ri.ri2py_dataframe(sample.rx2('mineral_sat')).apply(pd.to_numeric, errors='coerce')})
                # replace sample mineral_sat entry with None if there is no mineral saturation data.
                if(len(dict_sample_data['mineral_sat'].index) == 1 and dict_sample_data['mineral_sat'].index[0] == 'None'):
                    dict_sample_data['mineral_sat'] = None

            if get_redox:
                dict_sample_data.update(
                    {"redox": pandas2ri.ri2py_dataframe(sample.rx2('redox')).apply(pd.to_numeric, errors='coerce')})

            if get_charge_balance:
                dict_sample_data.update({"charge_balance": df_charge_balance.loc[sample.rx2('name')[0], :]})
            
            if get_ion_activity_ratios:
                try:
                    dict_sample_data.update(
                        {"ion_activity_ratios": pandas2ri.ri2py_dataframe(sample.rx2('ion_activity_ratios'))})
                except:
                    dict_sample_data['ion_activity_ratios'] = None
            
            if get_affinity_energy:
                dict_sample_data.update({"affinity_energy_raw": pandas2ri.ri2py_dataframe(
                    sample.rx2('affinity_energy_raw')).apply(pd.to_numeric, errors='coerce')})
                dict_sample_data.update(
                    {"affinity_energy": pandas2ri.ri2py_dataframe(sample.rx2('affinity_energy')).apply(pd.to_numeric, errors='coerce')})

            out_dict["sample_data"].update(
                {sample_data.names[i]: dict_sample_data})

        out_dict.update({"batch_3o": batch_3o})
        speciation = Speciation(out_dict)

        if report_filename != None:
            if ".csv" in report_filename[-4:]:
                out_dict["report"].to_csv(report_filename)
            else:
                out_dict["report"].to_csv(report_filename+".csv")

        if delete_generated_folders:
            self._delete_rxn_folders()
        
        if self.verbose > 0:
            print("Finished!")
        
        return speciation


    def create_data0(self,
                     filename,
                     filename_ss=None,
                     data0_formula_ox_name=None,
                     suppress_redox=[],
                     db="wrm",
                     exceed_Ttr=True,
                     grid_temps=[0.0100, 50.0000, 100.0000, 150.0000,
                                 200.0000, 250.0000, 300.0000, 350.0000],
                     grid_press="Psat",
                     infer_formula_ox=False,
                     generate_template=True,
                     template_name=None,
                     verbose=1):
        """
        Create a data0 file from a custom thermodynamic dataset.
        
        Parameters
        ----------
        filename : str
            Name of csv file containing thermodynamic data in the OBIGT format.
            
        filename_ss : str, optional
            Name of file containing solid solution parameters.
        
        data0_formula_ox_name : str, optional
            Name of supplementary file containing data0 parameters and inferred
            formula oxidation states. Ignored if `infer_formula_ox` is False.
            See `infer_formula_ox` for more detail.
        
        suppress_redox : list of str, default []
            Suppress equilibrium between oxidation states of listed elements
            (Cl, H, and O cannot be included).
        
        db : str, default "wrm"
            Desired three letter code of data0 output.

        exceed_Ttr : bool, default True
            Calculate Gibbs energies of mineral phases and other species
            beyond their transition temperatures?

        grid_temps : list of eight float, default [0.0100, 50.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000, 350.0000]
            Eight temperature values that make up the T-P grid.
        
        grid_press : list of float, default "Psat"
            Eight pressure values that make up the T-P grid. "Psat" for
            calculations along the liquid-vapor saturation curve.

        infer_formula_ox : bool, default False
            Create a supplementary file containing data0 parameters and
            inferred formula oxidation states? This option is useful for
            creating as many entries in the formula_ox column when creating a
            new supplementary file. Note that compounds like DySO4+ result in
            blank entries in formula_ox because the redox states of two
            elements, Dy and S, would have to be estimated together; S has many
            oxidation states and Dy's oxidation states are not hard-coded.
        
        generate_template : bool, default True
            Generate a CSV sample input template customized to this data0?
        
        template_name : str, optional
            Name of the sample input template file generated. If no name is
            supplied, defaults to 'sample_template_xyz.csv', where 'xyz' is
            the three letter code given to `db`. Ignored if `generate_template`
            is False.
        
        verbose : int, 0, 1, or 2, default 1
            Level determining how many messages are returned during a
            calculation. 2 for all messages, 1 for errors or warnings only,
            0 for silent.
        """
        
        # Check that thermodynamic database input files exist and are formatted
        # correctly.
        self._check_database_file(filename)
        if filename_ss != None:
            self.__file_exists(filename_ss)
        
        self.verbose = verbose
        
        if self.verbose >= 1:
            print("Creating data0.{}...".format(db), flush=True)
        
        if sum([T >= 10000 for T in grid_temps]):
            raise Exception("Grid temperatures must be below 10000 °C.")
        
        if isinstance(grid_press, list):
            if sum([T >= 10000 for T in grid_temps]):
                raise Exception("Grid pressures must be below 10000 bars.")
        
        template = pkg_resources.resource_string(
            __name__, 'data0.min').decode("utf-8")
        grid_temps = convert_to_RVector(grid_temps)
        grid_press = convert_to_RVector(grid_press)
        suppress_redox = convert_to_RVector(suppress_redox)
        
        if filename_ss == None:
            filename_ss = ro.r("NULL")
        if data0_formula_ox_name == None:
            data0_formula_ox_name = ro.r("NULL")
        if template_name == None:
            template_name = "sample_template_{}.csv".format(db)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r_create_data0 = pkg_resources.resource_string(
                __name__, 'create_data0.r').decode("utf-8")
            ro.r(r_create_data0)
            ro.r.main_create_data0(filename=filename,
                                   filename_ss=filename_ss,
                                   grid_temps=grid_temps,
                                   grid_press=grid_press,
                                   db=db,
                                   template=template,
                                   exceed_Ttr=exceed_Ttr,
                                   data0_formula_ox_name=data0_formula_ox_name,
                                   suppress_redox=suppress_redox,
                                   infer_formula_ox=infer_formula_ox,
                                   generate_template=generate_template,
                                   template_name=template_name,
                                   verbose=self.verbose)
    
        for warning in w:
            print(warning.message)
        
        if self.verbose > 0:
            print("Finished creating data0.{}.".format(db))