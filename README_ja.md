# LS-SPH-Fluid-Simulator

> **重要なお知らせ (2026-08-24):**
> 本リポジトリの更新は終了しました. 
> パフォーマンス向上のために，本プロジェクトはリファクタリングされ，新しいリポジトリに移行されました．
> 
> 👉 **新しいリポジトリはこちら:** [LS-SPH-Benchmarks](https://github.com/shobuzako-kensuke/LS-SPH-Benchmarks)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15709255.svg)](https://doi.org/10.5281/zenodo.15709255)

[English README](./README.md)

本リポジトリは，**Least Squares Smoothed Particle Hydrodynamics (LS-SPH) method** [[1](#ref1)] に基づく様々なシミュレーションコードを提供しています．  

以下の計算が可能です．

1. 流体計算のベンチマークテスト
   - Taylor-Green vortex
   - Lid-driven cavity flow
   - Boussinesq convection

2. 地球物理および工学分野における応用問題（現在開発中）  


## ⚙️ 動作環境

計算には **Fortran** を使用し，その可視化には **Python** を使用します．  
必要な環境は以下の通りです．

| カテゴリ | 必要な環境 | メモ |
|:---|:---|:---|
|OS |Unix系環境 (Linux, MacOS, WSLなど) | テスト環境：Windows Subsystem for Linux (WSL2)|
|コンパイラ| Intel Fortran| テスト環境：`ifx` (デフォルト)
|ビルド | `make` | Fortran ファイルのビルドに利用|
|可視化 | `Python` | テスト環境：`Python 3.12.0` (`matplotlib` 等の基本的なライブラリが必要)|
|動画作成| `ffmpeg` | Pythonファイルで必要|

> [!TIP]
> WSL2のインストール方法や，WSL2上でFortranやPythonを実行するための環境構築などには，拙筆のQiita記事をご参照ください．  
> [WSL2のインストールとアンインストール](https://qiita.com/zakoken/items/61141df6aeae9e3f8e36)  
> [WSL2によるgfortranとintel fortranの環境構築](https://qiita.com/zakoken/items/2a5e629020ce68f3efe1)  
> [WSL2によるPython3の環境構築](https://qiita.com/zakoken/items/8ddfda7267e7d95b3c46)  

> [!TIP]
> (Unix系環境の場合) `ffmpeg` がインストールされていない場合は, ターミナルを開いて `sudo apt install ffmpeg` を実行してください．


## 🖥️ 使い方

1. 該当ディレクトリに移動
2. ターミナルから `make` を実行し，Fortranファイルをコンパイル
3. 続けて，`./start_calculation` を実行し，計算を開始
4. 計算終了後, `python main.py` を実行し，動画等を作成

> [!NOTE]
> 各問題設定は，それぞれのディレクトリに置かれた **README.md** をご参照ください．


## 🧑‍💻 引用情報

本リポジトリを利用・参考にされる際は，以下の**2つの文献**を引用してください．

<a id="ref1">[1]</a>  
Shobuzako, K., Yoshida, S., Kawada, Y., Nakashima, R., Fujioka, S., & Asai, M. (2025).  
A generalized smoothed particle hydrodynamics method based on the moving least squares method and its discretization error estimation.  
*Results in Applied Mathematics*, 26, 100594. [https://doi.org/10.1016/j.rinam.2025.100594](https://doi.org/10.1016/j.rinam.2025.100594)  

[2]  
Shobuzako, K. (2025). *LS-SPH-Fluid-Simulator* (Version 1.1.0) [Computer software]. Zenodo.  
[https://doi.org/10.5281/zenodo.15709255](https://doi.org/10.5281/zenodo.15709255)


## 🤝 プロジェクトへの貢献

コードの改善，バグの報告，または新機能の追加等にご協力いただける方は，お気軽にプルリクエストを送信してください．


## 🪪 ライセンス

本プロジェクトは [MIT ライセンス](./LICENSE) に準拠しています．
