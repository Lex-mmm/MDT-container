from __future__ import annotations

import sys
import os
import argparse  # <-- new import

import logging
import time
import uuid
from decimal import Decimal
from datetime import datetime

from sdc11073.location import SdcLocation
from sdc11073.loghelper import basic_logging_setup
from sdc11073.mdib import ProviderMdib
from sdc11073.mdib import mdibbase
from sdc11073.provider import SdcProvider
from sdc11073.provider.components import SdcProviderComponents
from sdc11073.roles.product import ExtendedProduct
from sdc11073.wsdiscovery import WSDiscoverySingleAdapter
from sdc11073.xml_types import pm_qnames as pm
from sdc11073.xml_types import pm_types
from sdc11073.xml_types.dpws_types import ThisDeviceType
from sdc11073.xml_types.dpws_types import ThisModelType
from sdc11073 import certloader


# Define the network adapter for discovery
NETWORK_ADAPTER = "en0"
# ------------- Original SDC Provider Setup -------------
# The provider’s UUID is created from a base.
base_uuid = uuid.UUID('{cc013678-79f6-403c-998f-3cc0cc050230}')
my_uuid = uuid.uuid5(base_uuid, "12345")


class SDCProvider:
    def __init__(self, active=False):

        self.active = active
        if not self.active:
            pass
        else:



            self.CertFolder = os.path.join(os.getcwd(), "Comms/static/ssl")
            print(f"Certificate folder: {self.CertFolder}")
            #self.CertFolder = "Comms/ssl/"
            self.mdibLocation = os.path.join(os.getcwd(), "Comms/static/mdib_draeger.xml")

            parser = argparse.ArgumentParser(description="Start SDC Provider")
            parser.add_argument('--adapter', default='en0', help="Network adapter to use (default: en0)")
            args = parser.parse_args()
            NETWORK_ADAPTER = args.adapter

            # Start discovery on the specified network adapter
            my_discovery = WSDiscoverySingleAdapter(NETWORK_ADAPTER)
            my_discovery.start()
            
            try:
                with open("Comms/static/mdib_draeger.xml", "r") as draegerMdib:
                    self.mdibFrame = ProviderMdib.from_mdib_file(self.mdibLocation)
            except Exception as e:
                print(f"Error reading MDIB file: {e}")
                # Set the location context for discovery.
            

            self.my_location = SdcLocation(fac='UMCU', poc='ICU_3', bed='01')  # Facility, Department, Bed

            # Set model information for discovery.
            dpws_model = ThisModelType(
                manufacturer='Draeger',
                manufacturer_url='www.draeger.com',
                model_name='VirtualPatientTestDevice',
                model_number='1.0',
                model_url='www.draeger.com/model',
                presentation_url='www.draeger.com/model/presentation'
            )
            dpws_device = ThisDeviceType(
                friendly_name='TestDevice',
                firmware_version='Version1',
                serial_number='12345'
            )

            # Load SSL context.
            my_ssl_context = certloader.mk_ssl_contexts_from_folder(
                ca_folder = self.CertFolder,
                ssl_passwd='dummypass'
            )

            # Create the SDC provider device.
            specific_components = SdcProviderComponents(role_provider_class=ExtendedProduct)
            self.sdc_provider = SdcProvider(
                ws_discovery=my_discovery,
                epr=my_uuid,
                this_model=dpws_model,
                this_device=dpws_device,
                device_mdib_container=self.mdibFrame,
                specific_components=specific_components,
                ssl_context_container=my_ssl_context
            )



    def start(self):
        print("Starting SDC provider...")
        self.sdc_provider.start_all()
        self.sdc_provider.set_location(self.my_location)

        ### Start provider, and initialize metric updates
        all_metric_descrs = [c for c in self.mdibFrame.descriptions.objects if c.NODETYPE == pm.NumericMetricDescriptor]
        with self.mdibFrame.metric_state_transaction() as transaction_mgr:
            for metric_descr in all_metric_descrs:
                st = transaction_mgr.get_state(metric_descr.Handle)
                st.mk_metric_value()
                st.MetricValue.Value = Decimal(1.0)
                st.MetricValue.ActiveDeterminationPeriod = 1494554822450
                st.MetricValue.Validity = pm_types.MeasurementValidity.VALID
                st.ActivationState = pm_types.ComponentActivation.ON
    
    def stop(self):
        print("Stopping SDC provider...")
        pass

    def sendData(self, dataType, data, active=False):
        self.active = active
        if not self.active:
            pass
        else:
            if dataType == "vital_sign":
                # Send vital sign data to the SDC provider
                timestamp = data['time']
                for metric in data['values']:
                    print(f"Sending vital sign data: {metric} = {data['values'][metric]}")
                    if metric == "mon_sat":

                        try:
                            with self.mdibFrame.metric_state_transaction() as transaction_mgr:
                                spoO2_state = transaction_mgr.get_state("SpO2.Measuredvalue.2F.44A5")
                                spoO2_state.MetricValue.Value = Decimal(data['values'][metric])
                        except Exception as e:
                            print(f"Warning: Failed to update MDIB metric (SpO2): {e}")
        
            elif dataType == "alarm":
                # Send alarm data to the SDC provider
                timestamp = data['time']
                metric = data['parameter']
            if metric == "mon_sat":
                try:
                    with self.mdibFrame.alert_state_transaction() as transaction_mgr:
                        ac_state = transaction_mgr.get_state("Limit.DESAT.SpO2.Measuredvalue.2F.44A5.39663")
                        as_state = transaction_mgr.get_state("AS.Vis.NORMALPRIOCOLORLATCHED.Limit.DESAT.SpO2.Measuredvalue.2F.44A5.39663.11566")
                        ac_state.Presence = True
                        as_state.Presence = pm_types.AlertSignalPresence.ON
                except Exception as e:
                    print(f"Warning: Failed to update MDIB alert state: {e}")
            else:
                # Default case or additional conditions can be added here.
                pass