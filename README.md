# stl_operator
- stlデータからそれを囲む直方体を作成し,`csv`,`vtk`に変換する.
- 分割stl,一括stlに対応.
  - 分割stl : 末端部1つずつ分離した状態のstlのこと.
    - 実行した際,最初の質問で`y`と答えること
    - 読み込む/出力するファイルの名前は`stl_operator.f90`の`set_filename`で管理している.
  - 一括stl : 末端部をまとめて分離した状態のstl
    - 実行した際,最初の質問で`n`と答えること
- `Fortran`では`csv`の変換までを行っている.`vtk`の変換は`csv2vtk.ipynb`を使いましょう.(`python`フォルダにある)