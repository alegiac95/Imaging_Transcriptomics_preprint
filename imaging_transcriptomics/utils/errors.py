# Exceptions defined for the Virtual Histology project.
# The exceptions include most of the possible errors to be dispalyed by the input controls.


class InvalidFormatError(Exception):
    """Exception raised when the format of one of the input files in not correct.

    Attributes:
        errorFile -- file that is not in the correct format.
        message -- optional user overridden error message to display.
    """

    def __init__(self, errorFile, message="The provided file has an invalid format. Please use files in the .nii, .nii.gz format."):
        self.errorFile = errorFile
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The error was caused by the file {self.errorFile}."


class InvalidSizeError(Exception):
    """Exception raised when the size of the images is not correct.

    Attributes:
        errorFile -- file with the wrong size
        size -- size of the input image
        message -- optional user defined error message
    """

    def __init__(self, errorFile, size, message="The provided file has a wrong size."):
        self.errorFile = errorFile
        self.message = message
        self.size = size
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} The file {self.errorFile} has size: {self.size}"
