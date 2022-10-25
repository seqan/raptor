#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

# non-preprocessing version

label: index reads

baseCommand: [ raptor, build ]

inputs:
  sequences:
    type:
      type: array
      items:
        type: array
        items: File

  kmer_size:
    type: int?
    label: The k-mer size
    doc: |
      Default: 20. Value must be in range [1,32]
    inputBinding:
      prefix: --kmer

  window_size:
    type: int?
    label: The window size
    doc: |
      Default: k-mer size. Value must be a positive integer.
    inputBinding:
      prefix: --window

  index_size:
    type: string?
    label: The size in bytes of the resulting index.
    doc: |
      Default: 1k. Must be an integer followed by [k,m,g,t] (case insenstive)
    inputBinding:
      prefix: --size

  output_name:
    type: string
    label: name of the index file to produce
    default: raptor.hibf.index
    inputBinding:
      prefix: --output

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
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: input_bins_filepaths.txt
        entry: |
          ${
             var bins = "";
             for (var i = 0; i < inputs.sequences.length; i++) {
                var currentBin = inputs.sequences[i];
                for (var j = 0; j < currentBin.length; j++) {
                  bins += currentBin[j].path + " ";
                }
                bins += "\n";
             }
             return bins;
          }

arguments:
  - prefix: --threads
    valueFrom: $(runtime.cores)
  - input_bins_filepaths.txt

outputs:
  index:
    type: File    
    outputBinding:
      glob: $(inputs.output_name)
