
## Origin: University Medical Center Utrecht (UMCU), team SAS-ICU, the Netherlands
## Version: 1.0
## Status: Review completed for testing purposes - Operational use in a clinical setting is pending Quality Management approval (UMCU)

## First review: R. Zoodsma, UMCU, 2nd of April 2025
## Second review: L. van Loon, UMCU, 31st of March 2025
## Updated as of 2nd of Aptil 2025


## Description
# This script contains mathematical functions to be used in the Ascom model for the dynamic alarming module in the light of SAS-ICU.
# Functionalities are split dependent on a) data-availability, alongside b) mathematical derivations
# Note that this list is not exhaustive and can be updated as the model progresses


# Libraries
import pandas as pd
import datetime

# Generic assumptions of stated functions:
    # data = dictionary of timestamp with corresponding value in format data = {"timestamp": [timestamps], "value": [values]}
    # assume time series continuity, missing values are filled with NULL or NaN
    # assume data frequency 1hz (1/sec)


## Functions
# 1. DATA AVAILABILITY FUNCTIONALITIES
class AscomMathFunctions:
    def __init__(self):
        pass
        
    # 1.1 Function to remove missing values from a dataset
    def QUAL_remove_missing(self, data):
        if data:
            data = pd.DataFrame(data)
            data = data.dropna(subset=['values'])
            data.reset_index(drop=True, inplace=True)  # Reset the indices
            return data if len(data)>0 else None
        else:
            return None

    # 1.2 Function to check for quantitative data availability - percentage of data points in comparison to the framelength
    def QUAL_avail_quant(self, data, quantPerct=0.75): ## just for testing purposes: 75% availability. valuess may need adjustment dependent on framelength

        # determine the completeness
        framelength = (max(data["time"]) - min(data["time"])).seconds

        data = self.QUAL_remove_missing(data) ## remove missing values

        completeness = len(data["values"]) / framelength
        if completeness >= quantPerct:
            return True
        else: ## not enough datapoints in comparison to framelength
            print(f"Data quantity is insufficient: {completeness} of {quantPerct}")
            return False
        
    # 1.3 Function to check for data availability against the framelength, checking data spread in comparison to framelength
    def QUAL_avail_frame(self, data, spreadPerct=0.95): ## just for testing purposes: need 99% availability. Values may need adjustment dependent on framelength
        # determine length of the time sequence
        frameLength = (max(data["time"]) - min(data["time"])).seconds
        data = self.QUAL_remove_missing(data) ## remove missing values
        dataLength = (max(data["time"]) - min(data["time"])).seconds ## length now corrected for min/max timestamp without missing values
        # determine the completeness
        if frameLength == 0:
            #print("Frame length is 0, no data available")
            return False
        completeness = dataLength / frameLength
        if completeness >= spreadPerct:
            return True
        else:
            print(f"Data spread is insufficient: {completeness} of {spreadPerct}")
            return False
        

        

    # 2. MATHEMATICAL FUNCTIONALITIES
    # Overview of data quality checks embedded per function:
        # 1. Mean timeseries: quantitative data check
        # 2. Mean blood pressure: quantitative data check    --- EXEMPLARY FUNCTION
        # 3. Median: quantitative data check
        # 4. First derivative: quantitative data check, framelength coverage
        # 5. Second derivative: quantitative data check, framelength coverage
        # 6. Area-under-the-curve: quantitative data check, framelength coverage
        # 7. Integral horizontal: quantitative data check, framelength coverage
        # 8. Integral vertical: quantitative data check, framelength coverage
        # 9. Envelop: None, Quality checks are embedded in [6, 7, 8]

    # Functions:

    # 2.1 Function to calculate the mean of a list in a single parameter, over multiple timepoints
    # NB: In line with previous discussions, the 'mean' function is separately denoted 'timeseries'. 
    # This is to distinguish between the mean of a value list over multiple timestamps (this function) as opposed to the mean of several values at a single timestamp (e.g. mean blood pressure dependent on Systolic/Diastolic value).
    def MATH_mean_timeseries(self, data):

        if not self.QUAL_avail_quant(data):
            return None ## just quantitative, framelength coverage of less importance here.
        else: ## continue
            data = self.QUAL_remove_missing(data) ## remove missing values
            mean = sum(data["values"]) / len(data["values"])
            return mean

    # 2.2 Function to calculate the mean of blood pressure -> EXEMPLARY, just to show the difference with the timeseries mean function (2.1)
    def MATH_mean_bp(systolic, diastolic):

        if not systolic or not diastolic:
            return None
        else: ## continue
            mean = (systolic + 2*diastolic) / 3
            return mean

    # 2.3 Function to calculate the median of a list
    def MATH_median(self, data):

        if not self.QUAL_avail_quant(data):
            return None
        else: ## continue
            data = self.QUAL_remove_missing(data) ## remove missing values

            data["values"].sort()
            n = len(data["values"])
            if n % 2 == 0:
                median = (data["values"][n//2] + data["values"][n//2 - 1]) / 2
            else:
                median = data["values"][n//2]
            return median

    # 2.4 Function to calculate the first derivative of a list
    def MATH_derivative_first(self, data):

        if not self.QUAL_avail_frame(data) or not self.QUAL_avail_quant(data):
            return None ## not enough data
        
        else:
            data = self.QUAL_remove_missing(data) ## remove missing values
            derivative = [data["values"][i+1] - data["values"][i] for i in range(len(data["values"])-1)]
            return derivative

    # 2.5 Function to calculate the second derivative of a list
    def MATH_derivative_second(self, data):
        # ASSUMPTIONS:
        # data = list of numerical values
        # Assume data continuity: data is continuous and in sequence, missing values will yield empty derivative. 
        # Missing value strategies are separately to be handled

        ## CHECK DATA AVAILABILITY, 
        if not self.QUAL_avail_frame(data) or not self.QUAL_avail_quant(data):
            return None ## not enough data
        
        else: ## continue
            derivative = self.MATH_derivative_first(data) ## first derivative, missing values are already removed
            ## second derivative is the derivative of the first derivative
            sec_derivative = [derivative[i+1] - derivative[i] for i in range(len(derivative)-1)]
            return sec_derivative
        
    # 2.6 Function to calculate the area-under-the-curve of a value list
    def MATH_auc(self, data, referenceValMin, referenceValMax):
        # ADDITIONAL ASSUMPTIONS:
        # Missing values will result in the no additional coverage due to no impuration, thus effectively underestimating AUC. However, this is a separate issue.
        # As such, missing value strategies will need to be handled separately - dependent on further use of the metric. 

        # Reference value is the starting value above/below which AUC is calculated
        # Differentiate between a 'lower' and 'upper' reference value.
        if not self.QUAL_avail_frame(data) or not self.QUAL_avail_quant(data):
            return None ## not enough data
        
        else: ## continue
            ## calculate AUC with reference value with absolute values
            data = self.QUAL_remove_missing(data) ## remove missing values

            AUC_below = sum([abs(data["values"][i] - referenceValMin) for i in range(len(data["values"])) if data["values"][i] < referenceValMin])
            AUC_above = sum([abs(data["values"][i] - referenceValMax) for i in range(len(data["values"])) if data["values"][i] > referenceValMax])
            AUC_total = AUC_below + AUC_above
            return {"AUC_total": int(AUC_total), "AUC_below": int(AUC_below), "AUC_above": int(AUC_above)}

    # 2.7 Function for integral calculation of a value list: Horizontal (Lebesgue's integral)
    def MATH_integral_horizontal(self, data, referenceValMin, referenceValMax):
        # ADDITIONAL ASSUMPTIONS:
        # Missing values will result in the no additional coverage due to no impuration, thus effectively underestimating true integral value. However, this is a separate issue.
        # As such, missing value strategies will need to be handled separately - dependent on further use of the metric. 

        # Reference value is the starting value above/below which integral is calculated
        # Differentiate between a 'lower' and 'upper' reference value.

        if not self.QUAL_avail_frame(data) or not self.QUAL_avail_quant(data):
            return None ## not enough data
        
        else: ## continue
            data = self.QUAL_remove_missing(data) ## remove missing values
            # Check whether below or above reference range: Square the absolute value, then sum
            integral_below = sum( [abs(data["values"][i] - referenceValMin) **2 for i in range(len(data["values"])) if data["values"][i] < referenceValMin] ) 
            integral_above = sum( [abs(data["values"][i] - referenceValMax) **2 for i in range(len(data["values"])) if data["values"][i] > referenceValMax] ) 
            integral_total = integral_below + integral_above
            
            return {"integral_total": int(integral_total), "integral_below": int(integral_below), "integral_above": int(integral_above)}

    # 2.8 Function for integral calculation of a value list: Vertical (Riemann's integral)
    def MATH_integral_vertical(self, data, referenceValMin, referenceValMax):
        # ADDITIONAL ASSUMPTIONS:
        # Missing values will result in the no additional coverage due to no impuration, thus effectively underestimating true integral value. However, this is a separate issue.
        # As such, missing value strategies will need to be handled separately - dependent on further use of the metric. 

        # Reference value is the starting value above/below which integral is calculated
        # Differentiate between a 'lower' and 'upper' reference value.

        if not self.QUAL_avail_frame(data) or not self.QUAL_avail_quant(data):
            return None ## not enough data
        
        else: ## continue
            data = self.QUAL_remove_missing(data) ## remove missing values
            data = pd.DataFrame(data)
            # Relative values of time and value, relative to reference values
            data["relative_time"] = [abs(data["time"][i] - min(data["time"])).total_seconds() for i in range(len(data["time"]))]
            data["relative_value"] = [int(0)] * len(data["values"]) ## initialize relative value list
            i= 0 ##index
            for value in data["values"]:
                i += 1
                if value < referenceValMin:
                    data.loc[i, "relative_value"] = int(value - referenceValMin)
                elif value > referenceValMax:
                    data.loc[i, "relative_value"] = int(value - referenceValMax)
                else:
                    data.loc[i, "relative_value"] = int(0)
                    

            unique_values = range(int(min(data["relative_value"])), 
                                int(max(data["relative_value"])) +1 ) ## range of unique values
            # Initialize integral values
            integral_below, integral_above = 0, 0

            for unique_val in unique_values:
                if unique_val < 0:
                    timeLength = len([data["relative_time"][i] for i in range(len(data["values"])) if data["values"][i] <= unique_val])
                    integral_below += timeLength**2  # Square the integral
                elif unique_val > 0:
                    timeLength = len([data["relative_time"][i] for i in range(len(data["values"])) if data["values"][i] >= unique_val])
                    integral_above += timeLength**2  # Square the integral

            integral_total = integral_below + integral_above
            return {"integral_total": int(integral_total), "integral_below": int(integral_below), "integral_above": int(integral_above)}
            

    # 2.9 EXEMPLARY FUNCTION: Function for calculation of the 'envelope' of a parameter
    def MATH_envelop(self, data, referenceValMin, referenceValMax):
        ## Envelope = AUC + Riemann integral + Lebesgue integral
        ## Quality check is done in the parent functions
        AUC = self.MATH_auc(data, referenceValMin, referenceValMax)
        RIEMANN = self.MATH_integral_vertical(data, referenceValMin, referenceValMax)
        LEBESGUE = self.MATH_integral_horizontal(data, referenceValMin, referenceValMax)

        if AUC and RIEMANN and LEBESGUE:
            envelop_total = AUC["AUC_total"] + RIEMANN["integral_total"] + LEBESGUE["integral_total"]
            envelop_below = AUC["AUC_below"] + RIEMANN["integral_below"] + LEBESGUE["integral_below"]
            envelop_above = AUC["AUC_above"] + RIEMANN["integral_above"] + LEBESGUE["integral_above"]
            # Output of this function is to be checked against the envelop-threshold (based on clinical context) to derive a decision: Alarm or no alarm
            return {"envelop_total": int(envelop_total), "envelop_below": int(envelop_below), "envelop_above": int(envelop_above)}
        else:
            return None
