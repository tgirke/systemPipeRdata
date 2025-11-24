cwlVersion: v1.0
class: CommandLineTool
baseCommand: blastn
arguments:
  - prefix: -out
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName)_blast.txt
    position: 99
inputs:
  query:
    type: File
    inputBinding:
      prefix: -query
  database_name:
    type: File
    inputBinding:
      prefix: -db
  evalue:
    type: int
    inputBinding:
      prefix: -evalue
  outfmt:
    type: int
    inputBinding:
      prefix: -outfmt
  results_path:
    type: Directory
  SampleName:
    type: string
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName)_blast.txt
