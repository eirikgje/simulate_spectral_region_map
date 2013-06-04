program simulate_spectral_region_map
   use healpix_types
   use quiet_mapfile_mod
   use utils
   use rngmod
   use udgrade_nr
   implicit none

   character(len=512)   :: paramfile, inmap_fname
   character(len=2)     :: itext
   integer(i4b) :: num_components, nside, npix, num_regions, currpix
   integer(i4b) :: i, j, k, l, curr_region_num
   integer(i4b), allocatable, dimension(:, :)   :: pixel_regions
   real(dp), dimension(3)                    :: pixvec
   integer(i4b) :: unit, seed, numband, nlist, ordering, red_fac, ninv_nside
   integer(i4b) :: amp_nside
   integer(i4b), allocatable, dimension(:)      :: behaviour
   real(dp), allocatable, dimension(:)  :: ref_freq, param_init, param_disp, a2t
   real(dp), allocatable, dimension(:)  :: freq
   real(dp)     :: x, radius, sigma_0, mul_factor
   real(dp), allocatable, dimension(:, :)       :: outmap, param_map, tempmap
   real(dp), allocatable, dimension(:, :)       :: region_param, ampmap, invN
   integer(i4b), allocatable, dimension(:, :)       :: pixel_region
   integer(i4b), allocatable, dimension(:)      :: numpix, listpix

   real(dp), parameter          :: K_BOLTZMANN = 1.3806503d-23 !In SI-units
   real(dp), parameter          :: H_PLANCK = 6.626068d-34 !In SI-units
   real(dp), parameter          :: T_CMB = 2.72548 !In SI-units
   type(planck_rng)     :: rng_handle

   unit = getlun()
   call getarg(1, paramfile)
   
   call get_parameter(unit, paramfile, 'NSIDE', par_int=nside)
   call get_parameter(unit, paramfile, 'SEED', par_int=seed)
   call rand_init(rng_handle, seed)
   npix = 12 * nside ** 2
   call get_parameter(unit, paramfile, 'NUM_COMPONENTS', par_int=num_components)
   call get_parameter(unit, paramfile, 'NUM_REGIONS', par_int=num_regions)
   allocate(behaviour(num_components))
   allocate(ref_freq(num_components))
   allocate(param_init(num_components))
   allocate(param_disp(num_components))
   allocate(ampmap(0:npix - 1, num_components))
   do i = 1, num_components
      call int2string(i, itext)
      call get_parameter(unit, paramfile, 'BEHAVIOUR_' // itext, par_int=behaviour(i))
      call get_parameter(unit, paramfile, 'REF_FREQ_' // itext, par_dp=ref_freq(i))
      call get_parameter(unit, paramfile, 'CENTRAL_PARAM_' // itext, par_dp=param_init(i))
      call get_parameter(unit, paramfile, 'PARAM_DISPERSION_' // itext, par_dp=param_disp(i))
      call get_parameter(unit, paramfile, 'AMP_MAP_' // itext, par_string=inmap_fname)
      call get_parameter(unit, paramfile, 'AMP_MAP_NSIDE_' // itext, par_int=amp_nside)
      call get_parameter(unit, paramfile, 'AMP_MAP_MUL_FACTOR_' // itext, par_dp= mul_factor)
      call read_map(tempmap, ordering, trim(inmap_fname))
      if (ordering == 1) then
         call convert_ring2nest(amp_nside, tempmap(:, 1))
      end if
      if (amp_nside /= nside) then
         call udgrade_nest(tempmap(:, 1) * mul_factor, amp_nside, ampmap(:, i), nside)
      else
         ampmap(:, i) = tempmap(:, 1) * mul_factor
      end if
      deallocate(tempmap)
   end do
   ref_freq = ref_freq * 1d9
   call get_parameter(unit, paramfile, 'NUMBAND', par_int=numband)
   allocate(freq(numband))
   do i = 1, numband
      call int2string(i, itext)
      call get_parameter(unit, paramfile, 'FREQ_' // itext, par_dp=freq(i))
   end do
   freq = freq * 1d9
   allocate(a2t(numband))
   allocate(invN(0:npix - 1, numband))
   allocate(pixel_regions(0:npix-1, num_components))

   call get_parameter(unit, paramfile, 'NINV_NSIDE', par_int=ninv_nside)
   do i = 1, numband
      call int2string(i, itext)
      !Assuming WMAP-style noise files here
      call get_parameter(unit, paramfile, 'NINV_' // itext, par_string=inmap_fname)
      call get_parameter(unit, paramfile, 'SIGMA0_' // itext, par_dp=sigma_0)
      call read_map(tempmap, ordering, trim(inmap_fname))
      if (ordering == 1) then
         call convert_ring2nest(ninv_nside, tempmap(:, 2))
      end if
      if (ninv_nside > nside) then
         red_fac = (ninv_nside / nside) ** 2
         do j = 0, npix-1
            invN(j, i) = sigma_0 / red_fac * sqrt(sum(1 / tempmap(j * & 
               & red_fac:(j+1) * red_fac-1, 2)))
         end do
      else if (nside < ninv_nside) then 
         stop "This just won't do"
      else
         invN(:, i) = sigma_0 / sqrt(tempmap(:, 2))
      end if
      x = H_PLANCK * freq(i) / (K_BOLTZMANN * T_CMB)
      a2t(i) = (exp(x) - 1) ** 2 / (x ** 2 * exp(x))
      deallocate(tempmap)
   end do

   allocate(outmap(0:npix-1, numband))
   allocate(param_map(0:npix-1, num_components))
   allocate(pixel_region(0:npix-1, num_components))
   allocate(numpix(num_regions))
   allocate(listpix(0:npix-1))

   do i = 1, num_components
      numpix = 0
      do while (any(numpix == 0))
         pixel_region(:, i) = 1
         numpix(1) = npix
         do j = 2, num_regions
            curr_region_num = j
            !Find a random pixel that will be the start of the new region
            !Not sure if this includes all pixels - but rand_uni gives numbers
            !from 0 to 1, 0 and 1 not inclusive, so I think this is right.
            currpix = int(rand_uni(rng_handle) * npix)
            !We make a disc around this pixel.
            call pix2vec_nest(nside, currpix, pixvec)
            !If somehow we don't get any pixels in this region, do it again
!            numpix(j) = 0
            nlist = 0
            do while (nlist == 0)
               !pi/2 was chosen arbitrarily, can be bigger
               radius = rand_uni(rng_handle) * pi / 2
               call query_disc(nside, pixvec, radius, listpix, nlist, nest=1)
            end do
            do k = 0, nlist - 1
               if (numpix(pixel_region(listpix(k), i)) == 1) then
                  continue
               else
                  numpix(pixel_region(listpix(k), i)) = & 
                     & numpix(pixel_region(listpix(k), i)) - 1
                  pixel_region(listpix(k), i) = curr_region_num
                  numpix(curr_region_num) = numpix(curr_region_num) + 1
               end if
            end do
         end do
         print *, numpix
      end do
   end do

   allocate(region_param(num_components, num_regions))
   do i = 1, num_components
      do k = 1, num_regions
         region_param(i, k) = param_init(i) + rand_gauss(rng_handle) * param_disp(i)
      end do
   end do

   param_map = 0
   do i = 0, npix - 1
      do j = 1, numband
         do k = 1, num_components
            do l = 1, num_regions
               if (pixel_region(i, k) == l) then
                  param_map(i, k) = region_param(k, l)
               end if
            end do
         end do
      end do
   end do

   outmap = 0
   do i = 0, npix - 1
      do j = 1, numband
         do k = 1, num_components
            if (behaviour(k) == 1) then
               !Power law
               outmap(i, j) = outmap(i, j) + ampmap(i, k) * a2t(j) * & 
                  & (freq(j) / ref_freq(k)) ** param_map(i, k)
            end if
         end do
         outmap(i, j) = outmap(i, j) + rand_gauss(rng_handle) * invN(i, j)
      end do
   end do

   do i = 1, num_components
      call int2string(i, itext)
      call write_map(param_map(:, i), 2, 'parammaps_' // itext // '.fits')
      call write_map(pixel_region(:, i), 2, 'regionmaps_' // itext // '.fits')
   end do
   do i = 1, numband
      call int2string(i, itext)
      call write_map(outmap(:, i), 2, 'signalmaps_' // itext // '.fits')
      call write_map(1 / (invN(:, i) ** 2), 2, 'invN_' // itext // '.fits')
   end do
end program simulate_spectral_region_map
