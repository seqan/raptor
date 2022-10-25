#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: search reads

baseCommand: [ raptor, search ]

inputs:
  index:
    type: File
    inputBinding:
       prefix: --index
  queries:
    type: File
    inputBinding:
       prefix: --query
    format:
      - edam:format_1929  # FASTA
      - edam:format_1930  # FASTQ
  false_positive_rate:
    type: double?
    label: The false positive rate used for building the index.
    doc: |
      Default: 0.05. Value must be in range [0,1]
    inputBinding:
      prefix: --fpr
  error_tolerance:
    type: int?
    label: The number of errors to tolerate
    doc: |
      Default: 0. Value must be a positive integer or 0.
    inputBinding:
      prefix: --error
  output_name:
    type: string
    label: name of the search results to produce
    default: raptor.search.output
    inputBinding:
      prefix: --output
  p_max:
    type: double
    label: Used in the dynamic thresholding
    doc: |
      The higher p_max, the lower the threshold. Default: 0.15. Value must be in
      the range [0,1]
    inputBinding:
      prefix: --p_max

hints:
  SoftwareRequirement:
    packages:
      raptor:
        specs: [ https://bio.tools/raptor ]
        version: [ "2.0.0"]
  DockerRequirement:
    dockerPull: quay.io/biocontainers/raptor:2.0.0--h19e8d03_1

requirements:
  EnvVarRequirement:
    envDef:
      SHARG_NO_VERSION_CHECK: "1"

arguments:
  - prefix: --threads
    valueFrom: $(runtime.cores)

outputs:
  results:
    type: File
    label: Search results
    doc: | 
      Format: query id + tab + comma separated list of bin numbers that contain matches
    outputBinding:
      glob: $(inputs.output_name)

$namespaces:
  edam: http://edamontology.org/
