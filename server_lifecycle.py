#!/usr/bin/env python3
def on_server_loaded(server_context):
    ''' If present, this function is called when the server first starts. '''
    print("on_server_loaded")
    pass

def on_server_unloaded(server_context):
    ''' If present, this function is called when the server shuts down. '''
    print("on_server_unloaded")
    pass

def on_session_created(session_context):
    ''' If present, this function is called when a session is created. '''
    print("on_session_created")
    pass

def on_session_destroyed(session_context):
    ''' If present, this function is called when a session is closed. '''
    print("on_session_destroyed")
    pass
