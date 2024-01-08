# ASFVBiotypingTool
Note that this tool requires that BLAST and python are installed on your computer, and the custom zip file includes a list of genes and their related biotypes as well as a zip file (that must be unzipped) containing all files in the program. All files should be kept together when possible. The python file itself can be unpacked for the basic function for non-commandline users.

### Example Usage
This is an example of using the files located within this folder. The ASFVG file is an example one - the Representative_p72 file should always be used. **Make sure to unpack the zip file!**

    python .\ASFVBiotyping.py -q .\ASFVG.fa
