module stl_operator_m
    use,intrinsic :: iso_c_binding

    implicit none

    type coordinate_t
        double precision coordinate(3)
    end type

    type triangle_t
        type(coordinate_t) vertexID(3)
        double precision normal_vector(3)
    end type

    type content_t
        integer ID
        double precision value
    end type

    type group_of_norm_t
        integer, allocatable :: faceIDs(:)
        double precision max_coordinate(3), min_coordinate(3)
        double precision center(3), length(3), half_length(3)
    end type

    type merge_triangle_t
        type(coordinate_t) vertices(500) !決め打ち
    end type

    type re_merge_triangle_t
        type(coordinate_t) vertices(500)
    end type

    type append_t
        integer vertexID
    end type

    type re_append_t
        integer vertexID
    end type

    type box_t
        type(coordinate_t), allocatable :: node(:)
    end type

    type processor_t
    type(triangle_t),allocatable :: triangle(:)
    type(group_of_norm_t), allocatable :: group_norm(:)
    type(merge_triangle_t), allocatable :: merge_triangles(:)
    type(re_merge_triangle_t), allocatable :: re_merge_triangles(:)
    type(append_t), allocatable :: append_IDs(:)
    type(re_append_t), allocatable :: re_append_IDs(:)
    integer num_triangles, num_groups
    character(80) header
    end type

    type stl_t
        type(processor_t), allocatable :: processor(:)
        type(content_t), allocatable :: sorted_array(:)
        integer cnt_diff_norm
        double precision basic_vector(3)

        contains

        procedure read_binary_stl
        procedure output_stl_ascii
        procedure get_min_max_coordinate_of_group
        procedure output_box_list
        procedure output_coordinate_csv
        procedure make_triangle_groups
        procedure output_stl_triangle_group
        procedure output_csv_triangle_group
        procedure set_filename

    end type

    contains

    logical function is_open_binary_stl_file(filename) result(is_open)
        character(*), intent(in) :: filename
        integer istat, n_unit

        open(newunit = n_unit, file = filename, action  = "read",&
            status = "old", form = "unformatted", access  = "stream", iostat = istat)

            if(istat == 0) then
                is_open = .true.
            else
                is_open = .false.
            end if

        close(n_unit)

    end function

    subroutine read_filelist(path_boundary, array, array_size)
        use simpleFile_reader
        character(*), intent(in) :: path_boundary
        character(len=50), allocatable :: boundarynumber(:)
        character(len=50) :: file_number
        integer, allocatable :: array(:)
        integer :: line, number, status, cnt
        integer, intent(out) :: array_size

        call read_textRecord(path_boundary, boundarynumber)

        cnt = 1
        array_size = size(boundarynumber)

        allocate(array(array_size))

        do line = 1, size(boundarynumber)
            file_number = trim(boundarynumber(line))
            read(file_number, *, iostat=status) number
            if (status == 0) then
                array(cnt) = number
                cnt = cnt + 1
            end if
        end do

        array_size = cnt - 1
        if (array_size > 0) then
            array = array(1:array_size)
        else
            deallocate(array)
        end if

    end subroutine

    subroutine set_filename(self, filename_array, output_fname_stl, output_fname_txt)
        class(stl_t) self
        character(*), allocatable, intent(out) :: filename_array(:), output_fname_stl(:), output_fname_txt(:)
        character(len=1) :: input
        character(len=100) :: stl_name, stl_dirname, cond_dirname, input_filename, cond_fname, cond_path, num_str, index_str
        integer, allocatable :: num_array(:)
        integer :: array_size, index
        logical existance

            print*, 'DO you want to perform calculations on the separated STL file? (y/n)'
            read(5,*) input

            select case(input)
            case('y')
                print *, 'condition file name?'
                read(5,'(a)') cond_fname
                print *, 'directory name of condition file?'
                read(5,'(a)') cond_dirname
                print *, 'directory name of stls?'
                read(5,'(a)') stl_dirname

                cond_path = trim(trim(cond_dirname)//'/'//trim(cond_fname))

                call read_filelist(path_boundary = cond_path, array = num_array, array_size = array_size)

                allocate(filename_array(array_size))
                allocate(output_fname_stl(array_size))
                allocate(output_fname_txt(array_size))

                do index = 1, array_size
                    write(num_str, '(i0)') num_array(index)
                    input_filename = trim(trim(stl_dirname)//'/'//"inlet"//trim(adjustl(num_str))//".stl")
                    inquire(file=input_filename, exist=existance)
                    if(.not.existance) then
                        print*, 'File:[ ', input_filename, '] is not found.'
                        error stop
                    end if
                    filename_array(index) = input_filename
                    write(index_str, '(i0)') index
                    output_fname_stl(index) = trim("triangle_"//trim(adjustl(index_str))//".stl")
                    output_fname_txt(index) = trim("min_max_coord_"//trim(adjustl(index_str))//".txt")
                    ! print *, input_filename
                end do
                allocate(self%processor(array_size))

            case('n')
                print*, 'stl name?'
                read(5,'(A)') stl_name
                print*, 'Directory name of stl?'
                read(5,'(A)') stl_dirname
                input_filename = trim(trim(stl_dirname)//'/'//trim(stl_name))
                filename_array = [input_filename]
                inquire(file=input_filename, exist=existance)
                if(.not.existance) then
                    print*, 'File:[ ', input_filename, '] is not found.'
                    error stop
                end if
                allocate(self%processor(1))
                allocate(output_fname_txt(1))
                output_fname_txt(1) = "min_max_coord.txt"

            end select
    end subroutine

    subroutine append(array, new_value)
        real, allocatable :: array(:), hold_array(:)
        real new_value
        integer current_size

        current_size = size(array)
        allocate(hold_array(current_size))
        hold_array = array !一時的に値を保持

        deallocate(array)
        allocate(array(current_size + 1))
        array(1:current_size) = hold_array
        array(current_size + 1) = new_value
        deallocate(hold_array)
    end subroutine

    subroutine read_binary_stl(self, filename, fileID)
        class(stl_t) self
        character(*), intent(in) :: filename
        integer, intent(in) :: fileID
        integer n_unit, istat, faceID
        integer(c_int32_t) num_triangles
        real(c_float) normal_vector(3)        
        real(c_float) first_vertex(3),second_vertex(3),third_vertex(3)
        integer(c_int16_t) byte_count

        if(is_open_binary_stl_file(filename)) then
            print*, "Successfully open binary file !"
        else
            print*, "Error : Couldn't open binary file"
            stop
        end if

            open(newunit = n_unit, file = filename, action  = "read",&
                status = "old", form = "unformatted", access  = "stream", iostat = istat)

                read(n_unit) self%processor(fileID)%header, num_triangles
                self%processor(fileID)%num_triangles = int(num_triangles)
                allocate(self%processor(fileID)%triangle(self%processor(fileID)%num_triangles))

                do faceID = 1, self%processor(fileID)%num_triangles
                    read(n_unit, iostat = istat) &
                    normal_vector,first_vertex,second_vertex,third_vertex,byte_count
                    if (istat /= 0) exit
                    self%processor(fileID)%triangle(faceID)%normal_vector = dble(normal_vector)
                    self%processor(fileID)%triangle(faceID)%vertexID(1)%coordinate = dble(first_vertex)
                    self%processor(fileID)%triangle(faceID)%vertexID(2)%coordinate = dble(second_vertex)
                    self%processor(fileID)%triangle(faceID)%vertexID(3)%coordinate = dble(third_vertex)
                end do

            close(n_unit)

    end subroutine

    subroutine output_stl_ascii(self, filename,fileID)
        class(stl_t) self
        character(*), intent(in) :: filename
        integer, intent(in) :: fileID
        integer n_unit, faceID

        open(newunit = n_unit, file = filename, status = "replace")

            write(n_unit,"(a)") "solid test"

            do faceID = 1, self%processor(fileID)%num_triangles
                write(n_unit,'(4x,*(g0:," "))') "facet normal", self%processor(fileID)%triangle(faceID)%normal_vector
                write(n_unit,"(8x,a)") "outer loop"
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%processor(fileID)%triangle(faceID)%vertexID(1)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%processor(fileID)%triangle(faceID)%vertexID(2)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%processor(fileID)%triangle(faceID)%vertexID(3)%coordinate
                write(n_unit,"(8x,a)") "endloop"
                write(n_unit,"(4x,a)") "endfacet"
            end do

            write(n_unit,"(a)") "endsolid test"
            
        close(n_unit)

    end subroutine

    subroutine output_coordinate_csv(self, filename, fileID)
        class(stl_t) self
        character(*), intent(in) :: filename
        integer, intent(in) :: fileID
        integer n_unit, faceID

        open(newunit = n_unit, file = filename, status = "replace")
            do faceID = 1, self%processor(fileID)%num_triangles
                write(n_unit, '(*(g0:,"  "))') faceID
                write(n_unit, '(*(g0:," , "))') self%processor(fileID)%triangle(faceID)%vertexID(1)%coordinate
                write(n_unit, '(*(g0:," , "))') self%processor(fileID)%triangle(faceID)%vertexID(2)%coordinate
                write(n_unit, '(*(g0:," , "))') self%processor(fileID)%triangle(faceID)%vertexID(3)%coordinate
            end do
        close(n_unit)

    end subroutine

    subroutine output_stl_triangle_group(self, faceID, filename, fileID)
        class(stl_t), intent(in) :: self
        integer n_unit, faceID, vertexID
        character(*), intent(in) ::  filename
        integer, intent(in) :: fileID

        open(newunit = n_unit, file = "data/triangle/stl/"//trim(filename), status = "replace")
            write(n_unit,"(a)") "solid test"
            do vertexID = 1, self%processor(fileID)%append_IDs(faceID)%vertexID - 1, 3
                write(n_unit,'(4x,*(g0:," "))') "facet normal", self%processor(fileID)%triangle(faceID)%normal_vector
                write(n_unit,"(8x,a)") "outer loop"
                write(n_unit,'(12x,*(g0:," "))') "vertex", &
                self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", &
                self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID + 1)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", &
                self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID + 2)%coordinate
                write(n_unit,"(8x,a)") "endloop"
                write(n_unit,"(4x,a)") "endfacet"
            end do
        close(n_unit)
    end subroutine

    subroutine make_triangle_groups(self, fileID, output_fname)
        implicit none
        class(stl_t), intent(inout) :: self
        integer faceID, vertexID, mergevertexID, n_unit, appendfaceID, new_size, &
        appendvertexID_1, appendvertexID_2, i, mergevertexID_1, mergevertexID_2, mergefaceID, &
        test_faceID, test_vertexID, re_appendvertexID, re_mergevertexID, del_vertexID, dummyID,&
        update_appendID, update_vertexID, update_faceID, count, logical_test_ID
        real index_value, coordinate(3), merge_coordinate(3), test_coordinate(3), &
        re_appendcoordinate(3), replace_coordinate(3)
        logical, allocatable :: classification_not_completed(:)
        logical all_false
        character(100) filename
        integer, intent(in) :: fileID
        character(*), allocatable, intent(in) :: output_fname(:)

        allocate(self%processor(fileID)%merge_triangles(200)) !決め打ち．気管支末端部の面の数より多ければ良い
        allocate(self%processor(fileID)%append_IDs(200)) !上の数字と同じにする

        !初期配列値
        do faceID = 1, size(self%processor(fileID)%merge_triangles)
            do vertexID = 1, size(self%processor(fileID)%merge_triangles(faceID)%vertices)
                self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID)%coordinate = (/-99,-99,-99/) !ダミー値-99を代入
                self%processor(fileID)%merge_triangles(1)%vertices(1)%coordinate &
                = self%processor(fileID)%triangle(1)%vertexID(1)%coordinate
                self%processor(fileID)%merge_triangles(1)%vertices(2)%coordinate &
                = self%processor(fileID)%triangle(1)%vertexID(2)%coordinate
                self%processor(fileID)%merge_triangles(1)%vertices(3)%coordinate &
                = self%processor(fileID)%triangle(1)%vertexID(3)%coordinate
            end do
        end do

        do i = 1, size(self%processor(fileID)%append_IDs)
            self%processor(fileID)%append_IDs(i)%vertexID = 1
        end do

        appendfaceID = 2
        do faceID =1, self%processor(fileID)%num_triangles
            print*, "triangleID : ", faceID
            mergeloop:do mergefaceID = 1, appendfaceID - 1
                do mergevertexID = 1, size(self%processor(fileID)%merge_triangles(faceID)%vertices)
                    do vertexID = 1, 3
                        coordinate = self%processor(fileID)%triangle(faceID)%vertexID(vertexID)%coordinate
                        merge_coordinate &
                        = self%processor(fileID)%merge_triangles(mergefaceID)%vertices(mergevertexID)%coordinate
                        index_value = norm2(coordinate - merge_coordinate)
                        ! print*, index_value
                        if(index_value <= 0.6) then !一致する座標があるとき
                            ! print*, "same coordinate"
                            ! print*, "mergefaceID = ", mergefaceID
                            appendvertexID_1 = self%processor(fileID)%append_IDs(mergefaceID)%vertexID
                            do mergevertexID_1 = 1, 3
                                coordinate = self%processor(fileID)%triangle(faceID)%vertexID(mergevertexID_1)%coordinate
                                self%processor(fileID)%merge_triangles(mergefaceID)&
                                %vertices(appendvertexID_1)%coordinate &
                                = coordinate
                                self%processor(fileID)%append_IDs(mergefaceID)%vertexID &
                                = self%processor(fileID)%append_IDs(mergefaceID)%vertexID + 1
                                appendvertexID_1 = self%processor(fileID)%append_IDs(mergefaceID)%vertexID
                            end do
                            exit mergeloop
                        else if(vertexID == 3 .and. &
                                mergevertexID == size(self%processor(fileID)%merge_triangles(appendfaceID - 1)%vertices)&
                                .and. mergefaceID == appendfaceID -1) then
                            ! print*, "not same coordinate"
                            ! print*, "mergefaceID = ", appendfaceID
                            appendvertexID_2 = self%processor(fileID)%append_IDs(appendfaceID)%vertexID
                            do mergevertexID_2 = 1, 3
                                coordinate = self%processor(fileID)%triangle(faceID)%vertexID(mergevertexID_2)%coordinate
                                self%processor(fileID)%merge_triangles(appendfaceID)%vertices(appendvertexID_2)%coordinate &
                                = coordinate
                                self%processor(fileID)%append_IDs(appendfaceID)%vertexID &
                                = self%processor(fileID)%append_IDs(appendfaceID)%vertexID + 1
                                appendvertexID_2 = self%processor(fileID)%append_IDs(appendfaceID)%vertexID
                            end do
                            appendfaceID = appendfaceID + 1
                            exit mergeloop
                        end if
                    end do
                end do
            end do mergeloop
        end do

        !全て格納した後，かぶりがないか再度チェック
        allocate(self%processor(fileID)%re_merge_triangles(200)) !merge_trianglesと同じ数
        allocate(self%processor(fileID)%re_append_IDs(200)) !append_IDsと同じ数

        do faceID = 1, size(self%processor(fileID)%merge_triangles)
            do update_vertexID = 1, size(self%processor(fileID)%merge_triangles(faceID)%vertices)
                self%processor(fileID)%re_merge_triangles(faceID)%vertices(update_vertexID)%coordinate&
                = self%processor(fileID)%merge_triangles(faceID)%vertices(update_vertexID)%coordinate
            end do
        end do
        do update_appendID = 1, size(self%processor(fileID)%append_IDs)
            self%processor(fileID)%re_append_IDs(update_appendID)%vertexID &
            = self%processor(fileID)%append_IDs(update_appendID)%vertexID
        end do

        allocate(classification_not_completed(size(self%processor(fileID)%merge_triangles)))
        do faceID = 1, size(self%processor(fileID)%merge_triangles)
            if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 == 0) then
                classification_not_completed(faceID) = .false.
            else
                classification_not_completed(faceID) = .true.
            end if
        end do

        all_false = all(.not. classification_not_completed)

        do while(.not. all_false)
            do faceID = 2, size(self%processor(fileID)%merge_triangles)
                if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0) then
                    do test_faceID = 1, faceID - 1
                        do test_vertexID = 1, self%processor(fileID)%append_IDs(test_faceID)%vertexID - 1
                            do vertexID = 1, self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                                coordinate = &
                                self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID)%coordinate
                                test_coordinate &
                                = self%processor(fileID)%merge_triangles(test_faceID)%vertices(test_vertexID)%coordinate
                                index_value = norm2(coordinate - test_coordinate)
                                if(index_value == 0) then
                                    !検出されたvertexIDが3で割れるかどうかの場合分け
                                    if(mod(vertexID, 3) /= 0) then
                                        re_appendvertexID = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID
                                        do re_mergevertexID = 3*int(vertexID/3)+1, 3*int(vertexID/3 + 1)
                                            !一致groupにtriangleの座標を追加
                                            re_appendcoordinate &
                                            = self%processor(fileID)%merge_triangles(faceID)%vertices(re_mergevertexID)%coordinate
                                            self%processor(fileID)%re_merge_triangles(test_faceID)&
                                            %vertices(re_appendvertexID)%coordinate &
                                            = re_appendcoordinate
                                            self%processor(fileID)%re_append_IDs(test_faceID)%vertexID &
                                            = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID + 1
                                            re_appendvertexID = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID
                                        end do
                                        !元のgroupから追加したtriangleの座標を削除(置き換え)
                                        do del_vertexID = 3*int(vertexID/3)+1, &
                                        self%processor(fileID)%append_IDs(faceID)%vertexID - 4
                                            replace_coordinate &
                                            = self%processor(fileID)%merge_triangles(faceID)%vertices(del_vertexID + 3)%coordinate
                                            self%processor(fileID)%re_merge_triangles(faceID)&
                                            %vertices(del_vertexID)%coordinate = replace_coordinate
                                        end do
                                        !ダミー値を増やす&append_IDsの更新
                                        do dummyID = self%processor(fileID)%append_IDs(faceID)%vertexID - 3, &
                                        self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                                            self%processor(fileID)%re_merge_triangles(faceID)%vertices(dummyID)%coordinate &
                                            = (/-99,-99,-99/)
                                        end do
                                        self%processor(fileID)%re_append_IDs(faceID)%vertexID &
                                        = self%processor(fileID)%append_IDs(faceID)%vertexID - 3
                                    else
                                        re_appendvertexID = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID
                                        do re_mergevertexID = 3*int(vertexID/3)-2, 3*int(vertexID/3)
                                            !一致groupにtriangleの座標を追加
                                            re_appendcoordinate &
                                            = self%processor(fileID)%merge_triangles(faceID)%vertices(re_mergevertexID)%coordinate
                                            self%processor(fileID)%re_merge_triangles(test_faceID)&
                                            %vertices(re_appendvertexID)%coordinate &
                                            = re_appendcoordinate
                                            self%processor(fileID)%re_append_IDs(test_faceID)%vertexID &
                                            = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID + 1
                                            re_appendvertexID = self%processor(fileID)%re_append_IDs(test_faceID)%vertexID
                                        end do
                                        !元のgroupから追加したtriangleの座標を削除
                                        do del_vertexID = 3*int(vertexID/3)-2, &
                                        self%processor(fileID)%append_IDs(faceID)%vertexID - 4
                                            replace_coordinate &
                                            = self%processor(fileID)%merge_triangles(faceID)%vertices(del_vertexID + 3)%coordinate
                                            self%processor(fileID)%re_merge_triangles(faceID)%vertices(del_vertexID)%coordinate &
                                            = replace_coordinate
                                        end do
                                        !ダミー値を増やす&append_IDsの更新
                                        do dummyID = self%processor(fileID)%append_IDs(faceID)%vertexID - 3, &
                                        self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                                            self%processor(fileID)%re_merge_triangles(faceID)%vertices(dummyID)%coordinate &
                                            = (/-99,-99,-99/)
                                            self%processor(fileID)%re_append_IDs(faceID)%vertexID &
                                            = self%processor(fileID)%append_IDs(faceID)%vertexID - 3
                                        end do
                                    end if
                                    !元配列の更新
                                    do update_faceID = 1, size(self%processor(fileID)%merge_triangles)
                                        do update_vertexID = 1, size(self%processor(fileID)%merge_triangles(faceID)%vertices)
                                            self%processor(fileID)%merge_triangles(update_faceID)%&
                                            vertices(update_vertexID)%coordinate&
                                            = self%processor(fileID)%re_merge_triangles(update_faceID)%&
                                            vertices(update_vertexID)%coordinate
                                        end do
                                    end do
                                    do update_appendID = 1, size(self%processor(fileID)%append_IDs)
                                        self%processor(fileID)%append_IDs(update_appendID)%vertexID &
                                        = self%processor(fileID)%re_append_IDs(update_appendID)%vertexID
                                    end do
                                end if
                            end do
                        end do
                    end do
                end if
            end do
            do logical_test_ID = 1, size(self%processor(fileID)%merge_triangles)
                if(self%processor(fileID)%append_IDs(logical_test_ID)%vertexID - 1 > 10 .or. &
                self%processor(fileID)%append_IDs(logical_test_ID)%vertexID - 1 == 0) then
                    classification_not_completed(logical_test_ID) = .false.
                end if
            end do
            all_false = all(.not. classification_not_completed)
        end do

        count = 1
        if(size(self%processor) == 1) then
            do faceID = 1, size(self%processor(fileID)%merge_triangles)
                if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0)then
                        write(filename,'("triangle_",i0,".stl")') count
                    call output_stl_triangle_group(self, faceID, filename, fileID)
                    count = count + 1
                end if
                self%processor(fileID)%num_groups = count - 1
            end do
        else
            do faceID = 1, size(self%processor(fileID)%merge_triangles)
                if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0)then
                    call output_stl_triangle_group(self, faceID, output_fname(fileID), fileID)
                end if
            end do
            self%processor(fileID)%num_groups = 1
        end if
        
        print*, "num of face group = ", self%processor(fileID)%num_groups

        count = 1
        open(newunit = n_unit, file = "data/triangle_group_coordinate.csv", status = "replace")
            do faceID = 1, size(self%processor(fileID)%merge_triangles)
                if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0)then
                write(n_unit,'(*(g0:,"  "))') count
                    do vertexID = 1, self%processor(fileID)%append_IDs(faceID)%vertexID - 1 !ダミー値を除いて記述
                        write(n_unit, '(*(g0:," , "))') self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID)%coordinate
                    end do
                    count = count + 1
                end if
            end do
        close(n_unit)

        if(fileID == 1) then
            open(newunit = n_unit, file = "data/triange_group_num.txt", status = "replace")
                write(n_unit, '(A)') "faceID, num_vertices"
                do faceID = 1, size(self%processor(fileID)%merge_triangles)
                    write(n_unit, '(*(g0:," , "))') faceID, self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                end do
            close(n_unit)
        end if
    end subroutine

    subroutine output_csv_triangle_group(self, fileID)
        class(stl_t) self
        integer n_unit, faceID, vertexID, count
        character(100) filename
        integer, intent(in) :: fileID

        count = 1
        do faceID = 1, size(self%processor(fileID)%merge_triangles)
            if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0)then
                write(filename, '("triangle_",i0,".csv")') count
                open(newunit = n_unit, file = "data/triangle/csv/"//trim(filename))
                    do vertexID = 1, self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                        write(n_unit, '(*(g0:," , "))') self%processor(fileID)%merge_triangles(faceID)%vertices(vertexID)%coordinate
                    end do
                    count = count + 1
                close(n_unit)
            end if
        end do
    end subroutine

    subroutine get_min_max_coordinate_of_group(self, fileID, output_fname)
        class(stl_t) self
        type(coordinate_t) node(5)
        integer groupID, dotID, faceID, n_unit, vertexID, cnt, lim_vertexID
        double precision max_coord(3), min_coord(3)
        integer, intent(in) :: fileID
        character(*), allocatable, intent(in) :: output_fname(:)

        allocate(self%processor(fileID)%group_norm(self%processor(fileID)%num_groups))

        cnt = 1
        do faceID = 1, size(self%processor(fileID)%merge_triangles)
            if(self%processor(fileID)%append_IDs(faceID)%vertexID - 1 > 0) then
                lim_vertexID = self%processor(fileID)%append_IDs(faceID)%vertexID - 1
                self%processor(fileID)%group_norm(cnt)%max_coordinate(1) &
                = maxval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(1))
                self%processor(fileID)%group_norm(cnt)%max_coordinate(2) &
                = maxval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(2))
                self%processor(fileID)%group_norm(cnt)%max_coordinate(3) &
                = maxval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(3))
                self%processor(fileID)%group_norm(cnt)%min_coordinate(1) &
                = minval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(1))
                self%processor(fileID)%group_norm(cnt)%min_coordinate(2) &
                = minval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(2))
                self%processor(fileID)%group_norm(cnt)%min_coordinate(3) &
                = minval(self%processor(fileID)%merge_triangles(faceID)%vertices(1:lim_vertexID)%coordinate(3))

                max_coord = self%processor(fileID)%group_norm(cnt)%max_coordinate
                min_coord = self%processor(fileID)%group_norm(cnt)%min_coordinate
                self%processor(fileID)%group_norm(cnt)%center = (max_coord + min_coord)/2.d0
                self%processor(fileID)%group_norm(cnt)%length = max_coord - min_coord
                self%processor(fileID)%group_norm(cnt)%half_length = self%processor(fileID)%group_norm(cnt)%length /2.d0

                cnt = cnt + 1
            end if
        end do

        open(newunit = n_unit, file = "data/txt/"//trim(output_fname(fileID)), status = "replace")
            do groupID = 1, size(self%processor(fileID)%group_norm)
                write(n_unit,'(*(g0:," "))') "groupID = ", groupID
                write(n_unit,'(*(g0:," "))') self%processor(fileID)%group_norm(groupID)%max_coordinate
                write(n_unit,'(*(g0:," "))') self%processor(fileID)%group_norm(groupID)%min_coordinate
            end do
        close(n_unit)

    end subroutine

    subroutine output_box_list(self,array)
        class(stl_t) self
        integer n_unit, groupID, fileID, scale_factor
        character(*), allocatable, intent(in) :: array(:)

        !stlのスケールがmm設定の場合,csvではmで表示する
        scale_factor = 1000

        open(newunit = n_unit, file = "data/box_list.csv", status = "replace")
            write(n_unit,'(*(g0:," "))') "x,y,z,lx,ly,lz"
            do fileID = 1, size(array)
                do groupID = 1, size(self%processor(fileID)%group_norm)
                    write(n_unit,'(*(g0:,","))') self%processor(fileID)%group_norm(groupID)%center/1000, &
                    self%processor(fileID)%group_norm(groupID)%length/scale_factor
                end do
            end do
        close(n_unit)

    end subroutine

    subroutine get_min_max_coord(node, min_coord, max_coord)
        type(coordinate_t), intent(in) :: node(:)
        double precision, intent(out) :: min_coord(:), max_coord(:)

        min_coord(1) = minval(node(:)%coordinate(1))
        min_coord(2) = minval(node(:)%coordinate(2))
        min_coord(3) = minval(node(:)%coordinate(3))

        max_coord(1) = maxval(node(:)%coordinate(1))
        max_coord(2) = maxval(node(:)%coordinate(2))
        max_coord(3) = maxval(node(:)%coordinate(3))

    end subroutine

    recursive subroutine merge_sort(array, left_end, right_end)
        type(content_t), intent(inout) :: array(:)
        integer, intent(in) :: left_end, right_end
        type(content_t), allocatable :: work_array(:)
        integer mid, i,j,k

        allocate(work_array(size(array)))
        
        if(left_end < right_end) then
            mid = int((left_end + right_end)/2)
            call merge_sort(array, left_end, mid)
            call merge_sort(array, mid+1, right_end)
            
            do i = mid, left_end, -1
                work_array(i) = array(i)           
            end do

            do j = mid+1, right_end
                work_array(right_end-(j-(mid+1))) = array(j)
            end do

            i = left_end
            j = right_end

            do k = left_end, right_end
                if(work_array(i)%value < work_array(j)%value) then
                    array(k) = work_array(i)
                    i = i + 1
                else
                    array(k) = work_array(j)
                    j = j - 1
                end if
            end do

        end if

    end subroutine

    function append2list_int(list, element) result(after_list)
        integer, allocatable, intent(in) :: list(:)
        integer, intent(in) :: element
        integer, allocatable :: after_list(:)
        integer n

        if(.not.allocated(list)) then
            allocate(after_list(1))
            after_list(1) = element
        else
            n = size(list)
            allocate(after_list(n+1))
            after_list(:n) = list(:n)
            after_list(n+1) = element
        end if

    end function
        
end module