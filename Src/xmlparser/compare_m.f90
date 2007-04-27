module compare_m
use flib_dom
use compare_tol_m, only: findTol, sp, dp
use m_strings

implicit none

logical :: test_debug = .false.

public compare
private

!This is crude ....
integer ::max_char_len = 100
contains

!-------------------------------------------------------------------------

recursive subroutine compare(reference,modified)
use flib_dom
! Given two nodes lists this subroutine compares all the nodes, one by one.
type(fnode), pointer :: reference
type(fnode), pointer :: modified

!Internal vars.
logical                  :: has_children_ref
integer                  :: i,n_children_ref,error
type(fnodelist),pointer  :: sub_list_ref
type(fnode),pointer      :: sub_reference, mod_node

has_children_ref = hasChildNodes(reference)

if (should_be_checked(reference)) then
   !print *,"should"
   mod_node => find_corresponding_mod_node(reference,modified,error)
   if (.not. associated(mod_node))then
      print*,"Error:No corresponding node in the file to be checked!"
      stop
   else
      call compare_node_list(reference,mod_node)
   endif

else

   sub_list_ref   => getchildNodes(reference)
   n_children_ref =  getLength(sub_list_ref)
   do i=0,n_children_ref-1
      sub_reference => item(sub_list_ref,i)      
      call compare(sub_reference,modified) 
   enddo

endif

end subroutine compare

!--------------------------------------------------------------------------

function should_be_checked(ref) result(should)
type (fnode), pointer :: ref
logical :: should 


character(len=400)       :: tag,attr_name
integer                  :: i,n_kind,length_attr
type(string)             :: value
type(fnamedNodeMap),pointer :: attr_ref
type(fnode)         ,pointer :: attr

should = .false.

!Look for the properties of this node.
if(test_debug) print *,"------------------------------------------"

tag = getTagName(ref)

if(test_debug) print *,"Tag: ", "|",trim(tag),"|"

if (tag == "") then
   should = .false.
   return
endif



if (trim(tag) == "cml") then
   should = .false.
   !print*," should: cml false"
   return
elseif(trim(tag) == "metadata")then
   should = .false.
   !print*," should: metadata false"
   return
elseif(trim(tag) == "parameterList") then
   !print*," should: paramslist true"
   should = .true.
   return
endif

n_kind = getNodeType(ref)
if(test_debug) print *,"   Type:",n_kind
value = getNodeValue(ref)

if (len(trim(value)) > 0)  print *, "Value:","|",char(value),"|"

if(hasAttributes(ref) )then
   attr_ref => getAttributes(ref)
   length_attr = getlength(attr_ref)
   !if (length_attr > 0) print *," Has",length_attr," attr:" 
   do i=0,length_attr-1
      attr => item(attr_ref,i)
      attr_name=getNodeName(attr)
      if (len(attr_name) > 0) then
         !print*, "  attr name:", trim(attr_name)
         !print*, "  attr value:", trim(getNodeValue(attr))
         should = .true.
         return
      endif
   enddo
endif

end function should_be_checked

!--------------------------------------------------------------------------

recursive subroutine compare_node_list(reference,modified)
use flib_dom
! Given two nodes lists this subroutine compares all the nodes, one by one.
type(fnode), pointer :: reference
type(fnode), pointer :: modified

!Internal vars.
logical                  :: has_children_ref,has_children_mod
integer                  :: i,n_children_ref, n_children_mod
type(fnodelist),pointer  :: sub_list_ref, sub_list_mod 
type(fnode),pointer      :: sub_reference, sub_mod
type(string)             :: sn_ref,sn_mod


has_children_ref = hasChildNodes(reference)
has_children_mod = hasChildNodes(modified)

call compare_node(reference,modified)

if ( has_children_ref .and. has_children_mod) then
   sub_list_ref   => getchildNodes(reference)
   n_children_ref =  getLength(sub_list_ref)

   sub_list_mod   => getchildNodes(modified)
   n_children_mod =  getLength(sub_list_mod)
   
   if ( n_children_ref /= n_children_mod ) then
      !print *, " The nodes have a different number of children!"
      sn_ref = getNodeName(reference)
      sn_mod = getNodeName(modified)
      !print *, " Ref. node name=",char(sn_ref)
      !print *, " Mod. node name=",char(sn_mod)
      !stop
   else
      do i=0,n_children_ref-1
         sub_reference => item(sub_list_ref,i)
         sub_mod => item(sub_list_mod,i)
         call compare_node_list(sub_reference,sub_mod) 
      enddo
   endif
endif

end subroutine compare_node_list

!--------------------------------------------------------------------------

function find_corresponding_mod_node(ref,mod,error) result (mod_node)
! Given a reference node this subroutine finds the corresponding node
! in the modified list.
type(fnode), pointer     :: ref  !The reference node 
type(fnode), pointer     :: mod  !The list of modified nodes.
type(fnode), pointer     :: mod_node ! The modified node corresponding to the reference node.
integer, intent(out)     :: error !Error code.

!Internal vars.
character(len=400)       :: tag
integer                  :: i,n_elem
type(fnode), pointer     :: node
type(fnodelist),pointer  :: mod_list
logical                  :: corresponding

nullify(mod_node)

error = 0

tag = getTagName(ref)
if (tag /= "") then
   mod_list => getElementsByTagName(mod,tag)
else
   print *,"Error in find_corresponding ..."
   stop
endif
n_elem=getLength(mod_list)

if (n_elem == 1) then
   if(test_debug)print *,"find_corresponding: Only one node -> found!"
   mod_node => item(mod_list,0)
else
   !Loop over the elements of the list and compare the nodes
   !until we found the corresponding one.
   if(test_debug) print *, "Number of elements in mod list:",n_elem
   do i=0,n_elem-1
      node => item(mod_list,i)
      corresponding = check_node(ref,node)  
      if (corresponding) then
         mod_node => node
         exit
      endif
   enddo
endif

end function find_corresponding_mod_node

!--------------------------------------------------------------------------
 
recursive subroutine compare_node(ref,mod)
  use flib_dom

  type(fnode),pointer :: ref
  type(fnode),pointer :: mod

  integer :: i,n_type_ref, n_type_mod,length_ref,length_mod
  type(string) :: sn_ref,sn_mod,sv_ref, sv_mod
  character(len=50) ::  label,cn_parent_ref
  type(fnamedNodeMap),pointer :: attr_ref,attr_mod
  type(fnode), pointer :: c_node_ref, c_node_mod
  logical :: same_values = .true.
  

  !Check name
  sn_ref= getNodeName(ref)
  sn_mod= getNodeName(mod)

  !Metada nodes hold siesta version, compiler, timer info. Skipp the comparision.
  if (sn_ref == "metadata") return


  !if (test_debug) print *, "  Compare node: Ref,mod node names:",trim(char(sn_ref)),trim(char(sn_mod))

  if (char(sn_ref) /= "" .and. char(sn_mod) /= "") then
     sn_ref = clean_string(sn_ref)
     sn_mod = clean_string(sn_mod)
  endif

  if (char(sn_ref) /= char(sn_mod)) call dump_error(ref,mod,name=sn_ref)

  !Check node type
  n_type_ref = getNodeType(ref)
  n_type_mod = getNodeType(mod)

  !if (n_type_ref /= n_type_mod ) call dump_error(ref,mod,type=n_type_ref)

  !Check node value function should_be_checked. This part is the most complicate.
  sv_ref = getNodevalue(ref)
  sv_mod = getNodevalue(mod)

  sv_ref = clean_string(sv_ref)
  sv_mod = clean_string(sv_mod)

  !Check values
  !if (test_debug) print *," compare: values:",trim(char(sv_ref))
  
  if (len(sv_ref) > 0) then

     if (test_debug) then
        print *,"Name :","|",char(sn_ref),"|"
        print *,"Value:","|",char(sv_ref),"|", char(sv_mod),"|"
     endif

     !If the name is empty then look for parents name.
     if (sn_ref == "" .or. sn_ref=="#text")then
        call getParentNodeProperties(ref,cn_parent_ref)
        label = trim(cn_parent_ref)
     else
        label = trim(char(sn_ref))
     endif

     !Main comparision
     same_values = compare_values(ref,mod,label)

     if ( .not. same_values) then
        call dump_error(ref,mod,value=sv_ref)
     endif
  endif

  if (hasAttributes(ref) .and. hasAttributes(mod)) then

     attr_ref => getAttributes(ref)
     attr_mod => getAttributes(mod)

     length_ref = getlength(attr_ref)
     length_mod = getlength(attr_mod)

     if (test_debug) then
        print *, "   Ref, mod, attributes length:"
        print *, "     ",length_ref,length_mod
     endif

     if (length_ref /= length_mod) then
        call dump_error(ref,mod,attr=attr_ref)
     else
        do i=0,length_ref-1
           c_node_ref => item(attr_ref,i)
           c_node_mod => item(attr_mod,i)           
           call compare_node(c_node_ref, c_node_mod)
        enddo
     endif
  endif

!print *,"-----------------------"

endsubroutine compare_node

!-------------------------------------------------------------------------

recursive function check_node(ref,mod) result (same)
  use flib_dom

  type(fnode),pointer :: ref
  type(fnode),pointer :: mod
  logical :: same

  integer :: i,n_type_ref, n_type_mod,length_ref,length_mod
  type(string) :: sn_ref,sn_mod,sv_ref, sv_mod
  type(string) :: attr_name_ref, attr_value_ref, attr_value_mod, attr_name_mod
  type(fnamedNodeMap),pointer :: attr_ref,attr_mod
  type(fnode), pointer ::  c_node_ref, c_node_mod

  same = .true.

  !Check name
  sn_ref= getNodeName(ref)
  sn_mod= getNodeName(mod)

  if (test_debug)  print *, " Ref,mod node names:",char(sn_ref),"|",char(sn_mod)

  if (char(sn_ref) /= "" .and. char(sn_mod) /= "") then
     sn_ref = clean_string(sn_ref)
     sn_mod = clean_string(sn_mod)
  endif

  if (char(sn_ref) /= char(sn_mod)) then
     !print *," check nodes: different names"
     same = .false.
     return
  endif

  !Check node type
  n_type_ref = getNodeType(ref)
  n_type_mod = getNodeType(mod)

  if (n_type_ref /= n_type_mod ) then 
     !print *," check nodes: different types"
     same = .false.
     return
  endif

  !Compare values if they aren't numeric
  sv_ref = getNodeValue(ref)
  sv_mod = getNodeValue(mod)

  if ( .not. only_numbers(sv_ref) .and. .not. only_numbers(sv_mod))then
     if(test_debug) print *," check: values:",char(sv_ref),char(sv_mod)
     if (sv_ref /= sv_mod) then
        same = .false.
        return
     endif
  endif

  if (hasAttributes(ref) .and. hasAttributes(mod)) then

     attr_ref => getAttributes(ref)
     attr_mod => getAttributes(mod)

     length_ref = getlength(attr_ref)
     length_mod = getlength(attr_mod)

     if (test_debug) print *, "   Ref, mod, attributes length:",length_ref,length_mod

     if (length_ref /= length_mod) then
        same = .false.
        return
     endif

     if (length_ref > 0)then
        do i=0,length_ref-1
           c_node_ref => item(attr_ref,i)
           c_node_mod => item(attr_mod,i)           
          
           !Check attributes names and values.
           attr_name_ref = getNodeName(c_node_ref)
           attr_name_mod = getNodeName(c_node_mod)

           if(test_debug) print *," check: attr names:",trim(char(attr_name_ref)),"|",trim(char(attr_name_mod))

           if (attr_name_ref /= attr_name_mod) then
              same = .false.
              return
           endif
              
           attr_value_ref = getNodeValue(c_node_ref)
           attr_value_mod = getNodeValue(c_node_mod)

            if(test_debug) print *," check: attr values:",char(attr_value_ref),char(attr_value_mod)

           !if ( .not. only_numbers(attr_value_ref) .and. .not. only_numbers(attr_value_mod))then
              if (attr_value_ref /= attr_value_mod) then
                 same = .false.
                 return
              endif

           !endif

        enddo
     endif
  endif

if(test_debug)print *," check: same?",same
if(test_debug)print *,"-----------------------"

end function check_node
!----------------------------------------------

recursive subroutine getParentNodeproperties(node,name)
type(fnode), pointer                   :: node
character(len=50)                      :: name  !Name

!Internal vars
integer              :: i,length
type(fnode), pointer :: parent,attr
type(string)         :: n_attr
type(fnamedNodeMap),pointer :: attributes

parent => NULL()
parent => getParentNode(node)
if (associated(parent)) then
   name = getNodeName(parent) 
   !print*,"Begin ParentName: ",name

   !Not useful go up
   if (name == "" .or. name == "#text" .or. name == "scalar" .or. name == "array" &
        .or.name == "matrix")then ! .or. name == "property") then 
      call getParentNodeProperties(parent,name)
   elseif(name == "parameter" .or. name == "property")then

      attributes => NULL()
      if (associated(parent)) attributes => getAttributes(parent)

      if (.not. associated(attributes)) then
         call getParentNodeProperties(parent,name)
      else
         !Loop over the attributes and find the one with the name.
         length = getLength(attributes)
         do i=0,length-1
            attr => item(attributes,i)
            n_attr = getNodeName(attr)
            if ( n_attr /= "" .or. n_attr /= "#text" .or. n_attr /= "scalar" &
                 .or. n_attr == "title" .or. n_attr == "dictRef")then
               name = getNodeValue(attr)
               exit
            endif
         enddo
      endif
   endif
endif

end subroutine getParentNodeproperties

!-----------------------------------------------

subroutine dump_error(ref,mod,name,value,attr)
type(fnode), pointer :: ref,mod

type(string), intent(in), optional ::  name, value
!integer, intent(in), optional :: type
type(fnamedNodeMap),intent(in),optional :: attr

!Internal vars.
!integer :: n_type_ref, n_type_mod
type(string) :: sn_ref,sn_mod,sv_ref, sv_mod, snc_mod, snc_ref, svc_mod, svc_ref

logical ::  name_e = .false., value_e = .false., attr_e= .false.
character(len=50) :: n_parent_ref,n_parent_mod 
type(fnamedNodeMap),pointer :: attr_ref,attr_mod


!Find the type of error

!Find if the error is in the name:
if (present(name)) name_e = .true.

!Find if the error is in the value:
if (present(value)) value_e = .true.

!Find if the errror is in the attr.
if (present(attr)) attr_e = .true.


!Find all the info to make the differences report as
!complete as possible.

!Find parent properties.
call getParentNodeProperties(ref,n_parent_ref)
call getParentNodeProperties(mod,n_parent_mod)

!Name 
sn_ref = getNodeName(ref)
sn_mod = getNodeValue(mod)

snc_ref = clean_string(sn_ref)
snc_mod = clean_string(sn_mod)

!Type
!n_type_ref = getNodeType(ref)
!n_type_mod = getNodeType(mod)

!Values
sv_ref = getNodevalue(ref)
sv_mod = getNodevalue(mod)

svc_ref = clean_string(sv_ref)
svc_mod = clean_string(sv_mod)

!Attr.
attr_ref => getAttributes(ref)
attr_mod => getAttributes(mod)

!If the name is empty then use the parents name.
if (snc_ref == "" .or. snc_ref=="#text")then
   snc_ref = n_parent_ref
endif

if ( len(trim(sn_ref)) < 3 ) then
   call dump_error_heading(snc_ref,parent=n_parent_ref)
else 
   !if (str_ref == "scalar"
   call dump_error_heading(snc_ref)
endif

if (name_e) then  
   print *,"Different names:"
   print *,"  Ref: ",char(snc_ref)
   print *,"  Mod: ",char(snc_mod)
   stop
elseif(value_e)then
   print *, "Different values:"
   print *,"  Ref: ","|",char(svc_ref),"|"
   print *,"  Mod: ","|",char(svc_mod),"|"
   stop
elseif(attr_e)then
   print *, "Different attr:"
   !print *, char(sn_ref),char(sn_mod)
   stop
else
   
endif

stop "error"

end subroutine dump_error

!---------------------------------------------------

subroutine dump_error_heading (str,prop,parent)

type(string), intent(in) :: str
character(len=*), intent(in), optional :: prop,parent

print *,"There is an error in node: ","|",trim(char(str)),"|"
if (present(parent)) print *,"   which is a son of node: ", parent
end subroutine dump_error_heading

!---------------------------------------------------

function compare_values(ref,mod,label)
  ! Function wich compares two scalars (introduced in string format)
  ! If they do not differ then compare_scalar = .true.
  type(fnode), pointer            :: ref
  type(fnode), pointer            :: mod
  logical                         :: compare_values
  character(len=*), intent(in)    :: label


  !Internal vars.
  type(string) ::  n_ref, v_ref, n_mod, v_mod
  real(dp) :: tol

  n_ref = getNodeName(ref)
  v_ref = getNodeValue(ref)

  n_mod = getNodeName(mod)
  v_mod = getNodeValue(mod)

  v_ref = clean_string(v_ref)
  v_mod = clean_string(v_mod)

  tol = findTol(trim(label))

  !check if there are numbers in the string.
  if ( only_numbers(v_ref) .and. only_numbers(v_mod) )then 
     if (test_debug) print *, "v_ref,v_mod,label,tol",char(v_ref),char(v_mod),"|",label,"|",tol
     compare_values = compare_only_numbers(v_ref,v_mod,tol)
  else
     compare_values = compare_alpha(v_ref,v_mod)
  endif
 
end function compare_values

!---------------------------------------------------

function clean_string(str) result(cstr)
  use m_strings
  type(string),intent(in) :: str
  
  type(string) ::cstr
  
  !Remove the new lines

  cstr = str
  
  if (new_lines(cstr)) then
     if (test_debug) print *,"New lines!"
     call remove_new_lines(cstr)
  endif
  
  cstr = adjustl(cstr)
  cstr = trim(cstr)

end function clean_string

!---------------------------------------------------

function only_numbers(str)
  use m_strings
  type(string),intent(in)::str
  logical :: only_numbers

  character(len=53)::letters
  type(string) :: c_str
  integer :: length,position

  !Remove leading/trailing white spaces
  !print *,"original:","|",char(str),"|"
  
  letters="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ*"
  !Copy, we don't want to modify the xml file.
  c_str = str

  c_str = clean_string(c_str)

  length = len_trim(c_str)
  if (test_debug) print *,"|",char(c_str),"|"
  if (test_debug) print *,"length:",length
  
  position = scan(char(c_str),letters)

  if (position == 0)then
     only_numbers = .true.
  else
     only_numbers = .false.
  endif

   if (test_debug) print *,"only_numbers:",only_numbers
 
end function only_numbers
!---------------------------------------------------

function new_lines(str)
  use m_strings
  type(string),intent(in)::str
  logical :: new_lines
  
  integer :: position
  character :: nl
  
  nl = achar(10)

  if (len(str) > 0)then
     position = index(str,nl)
  else
     position = 0
  endif

  if (position == 0) then
     new_lines = .false.
  else
     new_lines = .true.
  endif

end function new_lines

!---------------------------------------------------
recursive subroutine remove_new_lines(str)
  use m_strings
  type(string),intent(inout)::str

  integer :: position
  character(len=1) :: nl
  
  nl = achar(10)
 
  do
     position = index(str,nl)
     if (position == 0)then
        exit
     else
        str = remove(str,position,position+1)
        call remove_new_lines(str)
     endif
  enddo

end subroutine remove_new_lines

!---------------------------------------------------
function compare_only_numbers(c_ref,c_mod,tol)
  type(string), intent(in) :: c_ref,c_mod
  real(dp), intent(in)         :: tol
  logical compare_only_numbers
  
  character :: ws
  integer   :: position
  real(dp)     :: r_ref,r_mod 
  type(string) :: ref,mod


  !Remove the new lines
  ref = c_ref
  mod = c_mod

  ref = clean_string(ref)
  mod = clean_string(mod)

  !Check if there are whitespaces:
  ws = achar(32)
  position = index(ref,ws)
  
  if (position == 0) then !There are no white spaces
     r_ref = str_to_scalar(ref)
     r_mod = str_to_scalar(mod)
     compare_only_numbers = compare_numbers(r_ref,r_mod,tol)
  else
     !There are whitespaces: compare array
     compare_only_numbers = compare_array(ref,mod,tol)     
  endif
end function compare_only_numbers

!---------------------------------------------------

function compare_alpha(c_ref,c_mod)
  type(string), intent(in) :: c_ref,c_mod
  logical compare_alpha
  
  
  !TODO: maybe something more elaborate?

  if (char(c_ref) /= char(c_mod)) then

     !print *,"  Alphanumeric strings don't match"
     !print *,"     Reference:",char(c_ref)
     !print *,"     Modified:",char(c_mod)
     !print *, "ref:",char(c_ref)
     !print *, "mod:",char(c_mod)
     compare_alpha = .true.
  else
     compare_alpha = .true.
  endif

end function compare_alpha

!---------------------------------------------------

function str_to_scalar(str) result(scal)
type(string),intent(in) :: str
real(dp) :: scal
integer  :: int

character(len=100) :: c_str

c_str = char(str)
!print*, "string:","|",char(str),"|"
if (is_integer(str)) then
   if (test_debug) print *,"string is integer:",c_str
   read(c_str,"(i3)"),int
   scal = real(int)
else
   read(c_str,"(f16.6)"),scal
endif
if (test_debug) print *,"Scalar:|",scal,"|"
end function str_to_scalar

!---------------------------------------------------

function is_integer(str)
  type(string), intent(in) :: str
  logical is_integer

  integer :: position
  character(len=2) :: po
  po = ".*"
  position = scan(str,po)
  if (position == 0)then
     is_integer = .true.
  else
     is_integer = .false.
  endif

end function is_integer

!---------------------------------------------------

function compare_numbers(ref,mod,tol)
  real(dp), intent(in) :: ref,mod,tol
  logical compare_numbers


  if (abs(abs(ref) - abs(mod)) .gt. tol)then
     compare_numbers = .false.
  else
     compare_numbers = .true.
  endif

end function compare_numbers
!---------------------------------------------------

function compare_array(ref,mod,tol)
  type(string), intent(in) :: ref, mod
  real(dp) , intent(in) :: tol
  logical :: compare_array

  integer :: i, length
  real(dp),allocatable,dimension(:)    :: a_ref,a_mod
  logical                          :: error

  call find_array(ref,mod, a_ref,a_mod,error)

  if (error) then
     compare_array = .false.
     return
  endif

  length = size(a_ref)

  compare_array = .true.

  if (length > 0)then
     do i=1,length
        
        if (abs(abs(a_ref(i))-abs(a_mod(i))) > tol) then   
           print *, "Compare array: error in element",i
           print *, trim(ref), " |", trim(mod)
           print *,"Values: ",abs(a_ref(i)),abs(a_mod(i))
           print *," Diff,tol: ", abs(a_ref(i))-abs(a_mod(i)),tol
           compare_array = .false.
           exit
        endif
     enddo
  endif

end function compare_array

!---------------------------------------------------

subroutine find_array(ref,mod,a_ref,a_mod,error)
  type(string), intent(in) :: ref, mod
  real(dp),allocatable,dimension(:) :: a_ref,a_mod
  logical, intent(out)          :: error

  integer      :: length = 1
  type(string) :: l_ref,l_mod,sub_ref,sub_mod
  integer :: pos_ref, pos_mod
  character :: ws
  character(len=2) ::ws2
  real(dp)         :: r_ref,r_mod

  error = .false.

  l_ref = ref
  l_mod = mod
  
  ws = achar(32)
  pos_ref = index(l_ref,ws)
  pos_mod = index(l_mod,ws)

  length = 0

  ws2="  "
  
  do
     pos_ref = index(l_ref,ws2)
     !print *,"pos_ref:",pos_ref
     if (pos_ref == 0)then
        exit
     else
        l_ref = remove(l_ref,pos_ref,pos_ref)
     endif
  enddo
 
  do
     pos_mod = index(l_mod,ws2)
     if (pos_mod == 0)then
        exit
     else
        l_mod = remove(l_mod,pos_mod,pos_mod)
     endif
  enddo
  
  pos_ref = index(l_ref,ws)
  pos_mod = index(l_mod,ws)

  do

     if (pos_ref == 0) exit

     sub_ref = extract(l_ref,1,pos_ref+1)
     sub_mod = extract(l_mod,1,pos_mod+1)
     
     !Get the scalars from the string.
     r_ref = str_to_scalar(sub_ref)
     r_mod = str_to_Scalar(sub_mod)
     
     !if (test_debug) print*, "Strings:",char(sub_ref),char(sub_mod)
     !if (test_debug) print*, "New scalars:",r_ref,r_mod

     !Remove the previous string
     l_ref = remove(l_ref,1,pos_ref)
     l_mod = remove(l_mod,1,pos_mod)
     
     !Store the values
     length = length + 1

     !Resize the array
     call resize(a_ref,length)
     call resize(a_mod,length)
      
     !store the new values
     a_ref(length) = r_ref
     a_mod(length) = r_mod

     !Look for the next ws.
     pos_ref = index(l_ref,ws)
     pos_mod = index(l_mod,ws)

     if (pos_ref > 0)then
        sub_ref = extract(l_ref,1,pos_ref+1)
        sub_mod = extract(l_mod,1,pos_mod+1)
     endif

  enddo
  

  if (len(sub_ref) > 0)then
     !print *, "one missing"
     r_ref = str_to_scalar(sub_ref)
     r_mod = str_to_Scalar(sub_mod)
     length = length + 1
     call resize(a_ref,length)
     call resize(a_mod,length)
      
     !store the new values
     a_ref(length) = r_ref
     a_mod(length) = r_mod
  endif

end subroutine find_array

!---------------------------------------------------

subroutine resize(a,new_length) 

  real(dp), dimension(:), allocatable, intent(inout) :: a
  integer, intent(in)            :: new_length

  !Internal vars
  real(dp), dimension(:), allocatable :: old_a
  integer                         :: length,i
  
  if (.not. allocated (a)) then
     allocate(a(1:new_length))
     a = 0.0
  else
     length = size(a)
     if (length > 0) then
        allocate(old_a(1:length))
        old_a = a
        deallocate(a)
        allocate(a(1:new_length))
        a = 0.0
        do i=1,size(old_a)
           a(i) = old_a(i)
        enddo
        deallocate(old_a)
     endif
  endif

end subroutine resize
!---------------------------------------------------

end module compare_m
