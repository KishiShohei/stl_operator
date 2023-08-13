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

    type box_t
        type(coordinate_t), allocatable :: node(:)
    end type

    type stl_t
        type(triangle_t),allocatable :: triangle(:)
        type(content_t), allocatable :: sorted_array(:)
        type(group_of_norm_t), allocatable :: group_norm(:)
        integer num_triangles, cnt_diff_norm
        double precision basic_vector(3)
        character(80) header

        contains

        procedure read_binary_stl
        procedure output_stl_ascii
        procedure output_dot_with_basic_vector
        procedure sort_dot_with_basic_vector
        procedure get_num_of_different_norm_vector
        procedure get_group_of_same_norm_vector
        procedure get_min_max_coordinate_of_group
        procedure output_box_list
        procedure output_box_vtk

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

    subroutine read_binary_stl(self, filename)
        class(stl_t) self
        character(*), intent(in) :: filename
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

            read(n_unit) self%header, num_triangles
            self%num_triangles = int(num_triangles)
            allocate(self%triangle(self%num_triangles))

            do faceID = 1, self%num_triangles
                read(n_unit, iostat = istat) &
                normal_vector,first_vertex,second_vertex,third_vertex,byte_count
                if (istat /= 0) exit
                self%triangle(faceID)%normal_vector = dble(normal_vector)
                self%triangle(faceID)%vertexID(1)%coordinate = dble(first_vertex)
                self%triangle(faceID)%vertexID(2)%coordinate = dble(second_vertex)
                self%triangle(faceID)%vertexID(3)%coordinate = dble(third_vertex)
            end do

        close(n_unit)

    end subroutine

    subroutine output_stl_ascii(self, filename)
        class(stl_t) self
        character(*), intent(in) :: filename
        integer n_unit, faceID

        open(newunit = n_unit, file = filename, status = "replace")

            write(n_unit,"(a)") "solid test"

            do faceID = 1, self%num_triangles
                write(n_unit,'(4x,*(g0:," "))') "facet normal", self%triangle(faceID)%normal_vector
                write(n_unit,"(8x,a)") "outer loop"
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%triangle(faceID)%vertexID(1)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%triangle(faceID)%vertexID(2)%coordinate
                write(n_unit,'(12x,*(g0:," "))') "vertex", self%triangle(faceID)%vertexID(3)%coordinate
                write(n_unit,"(8x,a)") "endloop"
                write(n_unit,"(4x,a)") "endfacet"
            end do

            write(n_unit,"(a)") "endsolid test"
            
        close(n_unit)

    end subroutine

    subroutine output_dot_with_basic_vector(self)
        class(stl_t) self
        integer n_unit, faceID

        ! 面数が合わない場合はここを変える.(任意の数)
        self%basic_vector(:) = (/0.1,0.2,0.3/)

        open(newunit = n_unit, file = "data/dot_with_basic_vector.txt", status = "replace")
            do faceID = 1, self%num_triangles
                write(n_unit,'(*(g0:," "))') faceID, &
                dot_product(self%basic_vector,self%triangle(faceID)%normal_vector)
            end do            
        close(n_unit)

    end subroutine

    subroutine sort_dot_with_basic_vector(self)
        class(stl_t) self
        type(content_t), allocatable :: array(:)
        integer n_unit, faceID, dotID

        allocate(array(self%num_triangles))
        
        do faceID = 1, self%num_triangles
            array(faceID)%ID = faceID
            array(faceID)%value = dot_product(self%basic_vector,self%triangle(faceID)%normal_vector)
        end do

        call merge_sort(array,1,self%num_triangles)

        open(newunit = n_unit, file = "data/sort_dot_with_basic_vector.txt", status = "replace")
            do dotID = 1, self%num_triangles
                write(n_unit,'(*(g0:," "))') array(dotID)%ID, array(dotID)%value
            end do
        close(n_unit)

    end subroutine

    subroutine get_num_of_different_norm_vector(self)
        class(stl_t) self
        integer n_unit, dotID, baseID, cnt, count_different_norm_vector

        allocate(self%sorted_array(self%num_triangles))

        open(newunit = n_unit, file = "data/sort_dot_with_basic_vector.txt", status = "old")
            do dotID = 1, self%num_triangles
                read(n_unit,*) self%sorted_array(dotID)%ID, self%sorted_array(dotID)%value
            end do
        close(n_unit)

        baseID = 1
        cnt = 2
        count_different_norm_vector = 1
        do dotID = cnt, self%num_triangles
            if(abs(self%sorted_array(baseID)%value - self%sorted_array(dotID)%value) <= 1.d-6) then
                cnt = cnt + 1
            else
                baseID = cnt
                cnt = cnt + 1
                count_different_norm_vector = count_different_norm_vector + 1    
            end if
        end do

        self%cnt_diff_norm = count_different_norm_vector
        print*, "count_different_norm_vector = ", self%cnt_diff_norm

    end subroutine

    subroutine get_group_of_same_norm_vector(self)
        class(stl_t) self
        integer baseID, cnt, dotID, faceID, cnt_diff_norm, n_unit, groupID

        allocate(self%group_norm(self%cnt_diff_norm))
        allocate(self%group_norm(1)%faceIDs(1))

        baseID = 1
        self%group_norm(1)%faceIDs(1) = self%sorted_array(baseID)%ID
        cnt = 2
        cnt_diff_norm = 1
        do dotID = cnt, self%num_triangles
            if(abs(self%sorted_array(baseID)%value - self%sorted_array(dotID)%value) <= 1.d-6) then
                cnt = cnt + 1

                faceID = self%sorted_array(dotID)%ID
                self%group_norm(cnt_diff_norm)%faceIDs &
                = append2list_int(self%group_norm(cnt_diff_norm)%faceIDs,faceID)
            else
                baseID = cnt
                cnt = cnt + 1
                cnt_diff_norm = cnt_diff_norm + 1

                faceID = self%sorted_array(dotID)%ID
                self%group_norm(cnt_diff_norm)%faceIDs &
                = append2list_int(self%group_norm(cnt_diff_norm)%faceIDs,faceID)
            end if
        end do

        open(newunit = n_unit, file = "data/same_norm_group.txt", status = "replace")
            do groupID = 1, size(self%group_norm)
                write(n_unit, '(*(g0:," "))') self%group_norm(groupID)%faceIDs
            end do
        close(n_unit)

    end subroutine

    subroutine get_min_max_coordinate_of_group(self)
        class(stl_t) self
        type(coordinate_t) node(5)
        integer groupID, dotID, faceID, n_unit, baseID
        double precision max_coord(3), min_coord(3)

        do groupID = 1, size(self%group_norm)
            ! グループ内の適当な座標で最大,最小座標を初期化
            baseID = self%group_norm(groupID)%faceIDs(1)
            min_coord(:) = self%triangle(baseID)%vertexID(1)%coordinate
            max_coord(:) = self%triangle(baseID)%vertexID(1)%coordinate

            do dotID = 1, size(self%group_norm(groupID)%faceIDs)
                faceID = self%group_norm(groupID)%faceIDs(dotID)

                node(1)%coordinate = self%triangle(faceID)%vertexID(1)%coordinate
                node(2)%coordinate = self%triangle(faceID)%vertexID(2)%coordinate
                node(3)%coordinate = self%triangle(faceID)%vertexID(3)%coordinate
                node(4)%coordinate = min_coord
                node(5)%coordinate = max_coord

                call get_min_max_coord(node, min_coord, max_coord)
                
            end do

            self%group_norm(groupID)%max_coordinate = max_coord
            self%group_norm(groupID)%min_coordinate = min_coord
            self%group_norm(groupID)%center = (max_coord + min_coord)/2.d0
            self%group_norm(groupID)%length = max_coord - min_coord
            self%group_norm(groupID)%half_length = self%group_norm(groupID)%length /2.d0

        end do

        open(newunit = n_unit, file = "data/min_max_coord.txt", status = "replace")
            do groupID = 1, size(self%group_norm)
                write(n_unit,'(*(g0:," "))') "groupID = ", groupID
                write(n_unit,'(*(g0:," "))') self%group_norm(groupID)%max_coordinate
                write(n_unit,'(*(g0:," "))') self%group_norm(groupID)%min_coordinate
            end do
        close(n_unit)

    end subroutine

    subroutine output_box_list(self)
        class(stl_t) self
        integer n_unit, groupID

        open(newunit = n_unit, file = "data/box_list.txt", status = "replace")
            write(n_unit,'(*(g0:," "))') "x,y,z,lx,ly,lz"
            do groupID = 1, size(self%group_norm)
                write(n_unit,'(*(g0:," "))') self%group_norm(groupID)%center, &
                self%group_norm(groupID)%length
            end do
        close(n_unit)

    end subroutine

    subroutine output_box_vtk(self)
        class(stl_t) self
        type(box_t), allocatable :: box(:)
        double precision, parameter :: &
        mat(3,8) = reshape([-1.d0, -1.d0, -1.d0, &
                             1.d0, -1.d0, -1.d0, &
                             1.d0,  1.d0, -1.d0, &
                            -1.d0,  1.d0, -1.d0, &
                            -1.d0, -1.d0,  1.d0, &
                             1.d0, -1.d0,  1.d0, &
                             1.d0,  1.d0,  1.d0, &
                            -1.d0,  1.d0,  1.d0], shape(mat))
        integer, allocatable :: vtk_mat(:)
        integer boxID, nodeID, n_unit, i

        allocate(box(size(self%group_norm)))
        do boxID = 1, size(box)
            allocate(box(boxID)%node(8))
        end do

        print*, mat(:,2)
        
        do boxID = 1, size(box)

            ! 閾値が0の場合, 幅をもたせる
            if(self%group_norm(boxID)%half_length(1) <= 1.d-4) &
            self%group_norm(boxID)%half_length(1) = 1.d-4
            if(self%group_norm(boxID)%half_length(2) <= 1.d-4) &
            self%group_norm(boxID)%half_length(2) = 1.d-4
            if(self%group_norm(boxID)%half_length(3) <= 1.d-4) &
            self%group_norm(boxID)%half_length(3) = 1.d-4

            do nodeID = 1, 8
                box(boxID)%node(nodeID)%coordinate &
                = self%group_norm(boxID)%center &
                  + mat(:,nodeID) * self%group_norm(boxID)%half_length
            end do
        end do

        open(newunit = n_unit, file = "data/box.vtk", status = "replace")
            write(n_unit,'(*(g0:," "))') '# vtk DataFile Version 2.0'
            write(n_unit,'(*(g0:," "))') 'Header'
            write(n_unit,'(*(g0:," "))') 'ASCII'
            write(n_unit,'(*(g0:," "))') 'DATASET UNSTRUCTURED_GRID'
            write(n_unit,'(*(g0:," "))') 'POINTS',size(box)*8,'float'
            do boxID = 1, size(box)
                do nodeID = 1, 8
                    write(n_unit,'(*(g0:," "))') box(boxID)%node(nodeID)%coordinate
                end do
            end do
            write(n_unit,'(*(g0:," "))') 'CELLS', size(box), size(box)*9
            allocate(vtk_mat(8))
            vtk_mat(:) = (/(i, i = 0,7)/)
            do boxID = 1, size(box)
                write(n_unit,'(*(g0:," "))') "8", vtk_mat + 8 * (boxID-1)
            end do
            write(n_unit,'(*(g0:," "))') 'CELL_TYPES', size(box)
            do boxID = 1, size(box)
                write(n_unit,'(*(g0:," "))') "12"
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