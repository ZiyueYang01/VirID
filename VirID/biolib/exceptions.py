class VirIDException(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)

class VirIDExit(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)

class FileNotFound(VirIDException):
    def __init__(self, message=''):
        VirIDException.__init__(self, message)


class DirNotFound(VirIDException):
    def __init__(self, message=''):
        VirIDException.__init__(self, message)  
        
class ExternalException(VirIDException):
    def __init__(self, message=''):
        VirIDException.__init__(self, message)

