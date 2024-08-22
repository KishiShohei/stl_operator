program main
    use stl_operator_m    
    implicit none

    character(256) input_filename, output_filename, csv_filename
    type(stl_t) stl

    input_filename = "data/outlet.stl"
    output_filename = "data/outlet_ascii.stl"
    csv_filename = "data/coordinate.csv"

    call stl%read_binary_stl(trim(input_filename))
    call stl%output_stl_ascii(trim(output_filename))
    call stl%output_coordinate_csv(trim(csv_filename))
    call stl%make_triangle_groups()
    call stl%output_csv_triangle_group()
    ! call stl%output_dot_with_basic_vector()
    ! call stl%sort_dot_with_basic_vector()
    ! call stl%get_num_of_different_norm_vector()
    ! call stl%get_group_of_same_norm_vector()

    !なんかセグメンテーション違反した．とりあえずコメントアウト
    ! call stl%get_min_max_coordinate_of_group()
    ! call stl%output_box_list()
    ! call stl%output_box_vtk()

end program