# ADRnet Training 
Adverse Drug Reaction prediction with improved data processing and expanded dataset training.

## Key Updates

### 1. Data Processing Pipeline
- **Chem Descriptor Calculation**:  
  Generated using OpenBabel for molecular fingerprinting (ECFP2, MACCS, mordred, RDkit, Morgan)

### 2. Extended Dataset Training
- **Original Training Data**:
  - AEOLUS database
  - Liu et al. dataset
- **New Training Data**:  
  - SIDER 4.1 (Side Effect Resource)
  - TOXRIC database

## Reference
```bibtex
@article{li2023generalized,
  author = {Li, Haoxuan and Hu, Taojun and Xiong, Zetong and Zheng, Chunyuan and Feng, Fuli and He, Xiangnan and Zhou, Xiao-Hua},
  title = {ADRNet: A Generalized Collaborative Filtering Framework Combining Clinical and Non-Clinical Data for Adverse Drug Reaction Prediction},
  year = {2023},
  isbn = {9798400702419},
  publisher = {Association for Computing Machinery},
  address = {New York, NY, USA},
  url = {https://doi.org/10.1145/3604915.3608813},
  doi = {10.1145/3604915.3608813},
  booktitle = {Proceedings of the 17th ACM Conference on Recommender Systems},
  pages = {682–687},
  numpages = {6},
  keywords = {adverse drug reaction, drug-ADR prediction, sided effect},
  location = {Singapore, Singapore},
  series = {RecSys '23}
}
