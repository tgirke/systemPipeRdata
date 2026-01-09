cwlVersion: v1.0
class: CommandLineTool
label: ""
doc: |
    written by Le Zhang
        11/2025
hints:
  SoftwareRequirement:
    packages:
    - package: fastqc
      version: [ 0.12.1 ]
      
baseCommand: [ fastqc ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.outdir) ]

arguments:
  - prefix: 
    valueFrom: $(inputs.fq1)
    position: 1
  
  - prefix: 
    valueFrom: $(inputs.fq2)
    position: 1
  
  - prefix: "--outdir"
    valueFrom: $(inputs.outdir)
    position: 3
    
  - prefix: "--threads"
    valueFrom: $(inputs.threads)
    position: 4
    
inputs:
  fq1:
    type: File
  fq2:
    type: File
  outdir:
    type: Directory
  threads:
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory
outputs:
  out_1:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq1.basename)_fastqc.html
  out_2:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq1.basename)_fastqc.zip
  out_3:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq2.basename)_fastqc.html
  out_4:
    type: File
    outputBinding:
      glob: $(inputs.outdir.path)/$(inputs.fq2.basename)_fastqc.zip
