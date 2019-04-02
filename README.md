# SCOtags
Some custom scripts related to the article Debray et al. 2019

# PhantomSpikesRemover.py
This script automates the removal of unusual subtitution rates that cause phantom spikes in Phylogenetic Informativeness profiles

Options:<br/>
- Required arguments:<br/>
  * `-d STR, --datapoints STR`<br/>
                        Path to reformated datapoints obtained through the PhyDesign application<br/>
  * `-r STR, --ratefile STR`<br/>
                        Path to the rate file obtained through the PhyDesign application<br/>
  * `-i STR, --infiles STR`<br/>
                        Path to the directory of individual alignments to clean. Must ends with a '/'.<br/>

- Optional arguments:<br/>
  * `-o OUTDIR, --outdir OUTDIR`<br/>
                        Output file name. Must ends with a '/'. Default is current directory `./`<br/>
  
  * `-h, --help`
                        Show the help message and exit<br/>


Usage example:<br/>
`python PhantomSpikesRemover.py -d DataPoints.reformated.txt -r RateFile.txt -i ./infiles/ -o ./oufiles/`
