import logging
import time
import uuid
import argparse
import os
import netifaces
import ipaddress
from sdc11073.xml_types import pm_types, msg_types
from sdc11073.xml_types import pm_qnames as pm
from sdc11073.xml_types.actions import periodic_actions
from sdc11073.wsdiscovery import WSDiscovery
from sdc11073.definitions_sdc import SdcV1Definitions
from sdc11073.consumer import SdcConsumer
from sdc11073.mdib import ConsumerMdib
from sdc11073 import observableproperties
from sdc11073.loghelper import basic_logging_setup
from sdc11073 import certloader
# This example shows how to implement a very simple SDC Consumer (client)
# It will scan for SDC Providers and connect to on well known UUID

# The provider we connect to is known by its UUID
# The UUID is created from a base
baseUUID = uuid.UUID('{cc013678-79f6-403c-998f-3cc0cc050230}')
device_A_UUID = uuid.uuid5(baseUUID, "12345")

# callback function that will be called upon metric updates from the provider
def on_metric_update(metrics_by_handle: dict):
    # we get all changed handles as parameter, iterate over them and output
    print(f"Got update on: {list(metrics_by_handle.keys())}")
    for handle in metrics_by_handle:
        print(handle.items())
        print(f"Handle: {handle}")

def on_alarm_update(alert_by_handle: dict):
    # we get all changed handles as parameter, iterate over them and output
    print(f"Got update on: {list(alert_by_handle.keys())}")

def set_ensemble_context(mdib: ConsumerMdib, sdc_consumer: SdcConsumer) -> None:
    # calling operation on remote device 
    print("Trying to set ensemble context of device A")
    # first we get the container to the element in the MDIB
    ensemble_descriptor_container = mdib.descriptions.NODETYPE.get_one(pm.EnsembleContextDescriptor)
    # get the context of our provider(client)
    context_client = sdc_consumer.context_service_client
    # start with empty operation handle and try to find the one we need
    operation_handle = None
    # iterate over all matching handles (can be 0..n)
    for op_descr in mdib.descriptions.NODETYPE.get(pm.SetContextStateOperationDescriptor, []):
        if op_descr.OperationTarget == ensemble_descriptor_container.Handle:
            operation_handle = op_descr.Handle
    # now we should have an operation handle to work with
    # create a new ensemble context as parameter to this operation
    new_ensemble_context = context_client.mk_proposed_context_object(ensemble_descriptor_container.Handle)
    new_ensemble_context.ContextAssociation = pm_types.ContextAssociation.ASSOCIATED
    new_ensemble_context.Identification = [
        pm_types.InstanceIdentifier(root="1.2.3", extension_string="SupervisorSuperEnsemble")]
    # execute the remote operation (based on handle) with the newly created ensemble as parameter
    response = context_client.set_context_state(operation_handle, [new_ensemble_context])
    result: msg_types.OperationInvokedReportPart = response.result()
    if result.InvocationInfo.InvocationState not in (msg_types.InvocationState.FINISHED,
                                                     msg_types.InvocationState.FINISHED_MOD):
        print(f'set ensemble context state failed state = {result.InvocationInfo.InvocationState}, '
              f'error = {result.InvocationInfo.InvocationError}, msg = {result.InvocationInfo.InvocationErrorMessage}')
    else:
        print(f'set ensemble context was successful.')


# main entry, will start to scan for the known provider and connect
# runs forever and consumes metrics everafter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Start SDC Consumer with specified network interface")
    parser.add_argument('--interface', default='en0', help="Network interface to use (default: en0)")
    parser.add_argument('--ssl-passwd', default="dummypass", help="SSL password for cert loading (default: dummypass)")
    args = parser.parse_args()

    # Convert interface name to IP address if needed
    try:
        # Check if the provided interface is already an IP address
        ipaddress.IPv4Address(args.interface)
        network_ip = args.interface
    except ipaddress.AddressValueError:
        # Get IP address from network interface name
        iface_data = netifaces.ifaddresses(args.interface)
        network_ip = iface_data[netifaces.AF_INET][0]['addr']

    # Determine current file's directory and set the ssl folder path (one folder deeper in 'ssl')
    current_dir = os.path.dirname(os.path.abspath(__file__))
    ssl_folder = os.path.join(current_dir, 'static/ssl')
    print(f"Using SSL folder: {ssl_folder}")
    basic_logging_setup(level=logging.INFO)
    my_discovery = WSDiscovery(network_ip)
    my_discovery.start()
    
    # we want to search until we found one device with this client
    found_device = False
    while not found_device:
        print('searching for sdc providers')
        services = my_discovery.search_services(types=SdcV1Definitions.MedicalDeviceTypesFilter)
        for one_service in services:
            print("Got service: {}".format(one_service.epr))
            if one_service.epr == device_A_UUID.urn:
                print("Got a match: {}".format(one_service))
                # Create SSL context container using the dynamic ssl folder
                my_ssl_context_container = certloader.mk_ssl_contexts_from_folder(
                    ca_folder=ssl_folder,
                    ssl_passwd=args.ssl_passwd
                )
                my_client = SdcConsumer.from_wsd_service(one_service, ssl_context_container=my_ssl_context_container)
                my_client.start_all(not_subscribed_actions=periodic_actions)
                my_mdib = ConsumerMdib(my_client)
                my_mdib.init_mdib()
                observableproperties.bind(my_mdib, metrics_by_handle=on_metric_update)
                #observableproperties.bind(my_mdib, alert_by_handle=on_alarm_update)
                found_device = True
                set_ensemble_context(my_mdib, my_client)
    
    while True:
        time.sleep(1)