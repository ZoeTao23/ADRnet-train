# ADR-Net Enhanced Training Framework
Enhanced framework for Adverse Drug Reaction prediction with improved data processing and expanded dataset training.

## Key Updates

### 1. New Data Processing Pipeline
- **ECFP Descriptor Calculation**:  
  Generated using OpenBabel for molecular fingerprinting
- **CACTVS Key Matching**:  
  Aligned with PubChem compound databases

### 2. Extended Dataset Training
- **Original Training Data**:
  - AEOLUS database
  - Liu et al. dataset
- **New Training Data**:  
  - SIDER 4.1 (Side Effect Resource)
  - TORIC 

## Reference
```bibtex
@article{li2023generalized,
  title        = {A Generalized Collaborative Filtering Framework Combining Clinical and Non-Clinical Data for Adverse Drug Reaction Prediction},
  author       = {Li, Haoxuan and Hu, Taojun and Xiong, Zetong and Zheng, Chunyuan and Feng, Fuli and He, Xiangnan and Zhou, Xiao-Hua},
  journal      = {[Journal Name]},
  volume       = {[Volume]},
  number       = {[Issue]},
  pages        = {[Page Range]},
  year         = {2023},
  publisher    = {[Publisher]},
  note         = {Manuscript submitted for publication}
}
