program main
    use stl_operator_m    
    implicit none

    character(256) input_filename, output_filename
    type(stl_t) stl

    input_filename = "data/test.stl"
    output_filename = "data/test_ascii.stl"

    call stl%read_binary_stl(trim(input_filename))
    call stl%output_stl_ascii(trim(output_filename))
    call stl%output_dot_with_basic_vector()
    call stl%sort_dot_with_basic_vector()
    call stl%get_num_of_different_norm_vector()
    call stl%get_group_of_same_norm_vector()
    call stl%get_min_max_coordinate_of_group()
    call stl%output_box_list()
    call stl%output_box_vtk()

end program