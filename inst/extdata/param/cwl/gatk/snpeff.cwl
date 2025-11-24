################################################################
##                        snpEf.cwl                       ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
label: ""
doc: |
    written by Le Zhang
        11/2025
hints:
  SoftwareRequirement:
    packages:
    - package: snpeff
      version: [ 5.3 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ java ]

arguments:
  - prefix: 
    valueFrom: "-Xmx30g -Xms16g -Djava.io.tmpdir=/tmp"
    position: 1
  
  - prefix: -jar
    valueFrom: snpEff/snpEff.jar 
    position: 2 
  
  - prefix: ann
    valueFrom: ""
    position: 3
  
  - prefix: -noStats
    valueFrom: ""
    position: 4
    
  - prefix: -q
    valueFrom: ""
    position: 5
    
  - prefix: 
    valueFrom:  $(inputs.genome_name)
    position: 6
    
  - prefix: 
    valueFrom:  $(inputs.filtered_vcf)
    position: 7
    
################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  gatk_java_options:
    type: string

  SampleName:
    type: string

  results_path:
    type: Directory

  filtered_vcf:
    type: File
  
  genome_name:
    type: string

stdout:  $(inputs.results_path.basename)/$(inputs.SampleName)_ann.vcf

outputs:
  ann_vcf:
    type: stdout

