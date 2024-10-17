program main
    use stl_operator_m    
    implicit none

    character(256) input_filename, output_filename, csv_filename, Boundaryfname
    type(stl_t) stl
    character(50), allocatable :: filename_array(:), output_fname_stl(:), output_fname_txt(:)
    integer caseID

    call stl%set_filename(filename_array,output_fname_stl,output_fname_txt)

    do caseID = 1, size(filename_array)
        call stl%read_binary_stl(trim(filename_array(caseID)),caseID)
        ! call stl%output_stl_ascii(trim(output_filename),caseID)
        ! call stl%output_coordinate_csv(trim(csv_filename),caseID)
        call stl%make_triangle_groups(caseID, output_fname_stl)
        ! call stl%output_csv_triangle_group(caseID)

        call stl%get_min_max_coordinate_of_group(caseID, output_fname_txt)
    end do

    call stl%output_box_list(filename_array)

end program