# Reproducible Single-Cell Analysis

This repository documents our attempt to reproduce the results of the paper ["Dictionary Learning for Integrative, Multimodal, and Scalable Single-Cell Analysis"](https://doi.org/10.1038/s41587-023-01767-y) as part of a seminar. Despite extensive efforts, we encountered significant challenges in replicating the results.

## Challenges and Difficulties

### Package Version Conflicts
- **Problem**: Conflicts between package versions led to numerous issues with package installations and functionality.
- **Details**: Specific packages required for the analysis had version conflicts that hindered their proper installation and operation.

### Sessions Crashing Due to Memory Limit
- **Problem**: R sessions frequently crashed due to memory constraints.
- **Details**: Even after attempting various memory management strategies, the sessions were consistently killed due to insufficient memory.

### Reproduction Failures
- **Problem**: Following the vignette and attempting to reproduce the results failed consistently.
- **Details**: Sessions were terminated each time we attempted to run the code, primarily due to high memory usage.

### Empirical Findings
- **Finding**: A minimum of 64GB of RAM was empirically determined to be necessary to run the code effectively.

## Actions Taken

- **Downsampled the Data**: Attempted to reduce the size of the data to fit within available memory limits.
- **Divided Data into Batches**: Processed data in smaller batches to manage memory usage better.
- **Allocated More Memory to R Session**: Tried increasing memory allocation to the R session.
- **Used Google Colab**: Leveraged Google Colab with ~32 GB RAM, but still faced issues due to memory limits.
- **Progress**: Despite these efforts, the session was killed every time due to insufficient memory.

## Conclusion

- **Difficulty in Reconstructing Paper's Data**: It was challenging to replicate the figures and results from the paper.
- **Significant Computational Resource Demand**: The analysis required substantial computational resources. The minimum hardware specification identified was 64GB of RAM.
- **Conflicting Packages**: Package versions and dependencies were unclear and contributed to the difficulties.

## Advice to the Authors

- **Clear Documentation**: Provide comprehensive documentation outlining the system requirements and steps needed to reproduce the results.
- **Minimum System Requirements**: Clearly list the minimal system requirements. Based on our findings, a minimum of 64GB of RAM is necessary.
- **Reproducibility Guidelines**: Offer detailed reproducibility guidelines to help users replicate the results effectively.
- **Community Support**: Engage with the community to provide support and address issues that users may face.
- **Version Control**: Use version control for both code and dependencies to ensure consistency and reproducibility.

## Reference

For further information, please refer to the original paper: ["Dictionary Learning for Integrative, Multimodal, and Scalable Single-Cell Analysis"](https://doi.org/10.1038/s41587-023-01767-y).

## Acknowledgments

We appreciate the effort that went into the original research and hope these insights will help in improving the reproducibility of the results.

For further information or to contribute, please refer to the discussions and documentation within this repository.
