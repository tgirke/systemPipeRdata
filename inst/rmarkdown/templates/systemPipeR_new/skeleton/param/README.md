# param folder

  - Stores non-CWL parameter files such as: *.param, *.tmpl and *.run.sh. These files are only required for backwards compatibility to run old workflows using the previous custom command-line interface.
  - param/cwl/: This subdirectory stores all the CWL parameter files. To organize workflows, each can have its own subdirectory, where all CWL param and input.yml files need to be in the same subdirectory.

The following example copies all the param templates files into new project template:

```
file.copy(systemPipeRdata::pathList()$paramdir, to = ".", recursive=TRUE, overwrite = TRUE)
```