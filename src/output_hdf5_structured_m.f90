   ! ! Logarithmic radial mesh
   ! call hdf_create_group(file_id, 'log_mesh')
   ! call hdf_open_group(file_id, 'log_mesh', log_mesh_group_id)
   ! call hdf_write_attribute(log_mesh_group_id, 'a', 0.01_dp * log(rr(101) / rr(1)))
   ! call hdf_write_attribute(log_mesh_group_id, 'b', rr(1))
   ! call hdf_write_dataset(log_mesh_group_id, 'r', rr)
   ! call hdf_set_data_scale(log_mesh_group_id, 'r', 'r (a.u.)')
   ! call hdf_write_attribute(log_mesh_group_id, 'r', 'units', 'Bohr')
   ! call hdf_write_attribute(log_mesh_group_id, 'r', 'description', 'r(i) = b * exp(a * (i - 1))')

   ! ! All-electron atom
   ! call hdf_create_group(file_id, 'ae_atom')
   ! call hdf_open_group(file_id, 'ae_atom', group_id)
   ! call hdf_write_dataset(group_id, 'v_tot', vfull)
   ! call hdf_attach_data_scale(log_mesh_group_id, 'r', group_id, 'v_tot')
   ! f_tmp(:) = rr(:)**2 * rhotae
   ! call hdf_write_dataset(file_id, 'rho_val', f_tmp)
   ! call hdf_attach_data_scale(log_mesh_group_id, 'r', group_id, 'rho_val')
   ! f_tmp(:) = rr(:)**2 * rhoc
   ! call hdf_write_dataset(file_id, 'rho_core', f_tmp)
   ! call hdf_attach_data_scale(log_mesh_group_id, 'r', group_id, 'rho_core')
   ! do l1 = 1, lmax + 1
   !    ll = l1 - 1

   ! end do
   ! call hdf_close_group(group_id)

   ! ! Pseudoatom
   ! call hdf_create_group(file_id, 'pseudo_atom')
   ! call hdf_open_group(file_id, 'pseudo_atom', group_id)
   ! call hdf_write_dataset()
   ! f_tmp(:) = rr(:)**2 * rho
   ! call hdf_write_dataset(group_id, 'rho_val', f_tmp)
   ! call hdf_attach_data_scale(log_mesh_group_id, 'r', group_id, 'rho_val')


   ! ! Pseudopotentials
   ! call hdf_create_group(file_id, 'pseudopotentials')
   ! call hdf_open_group(file_id, 'pseudopotentials', group_id)
   ! call hdf_write_dataset(group_id, 'v_sl', vp(:, 1:lmax + 1))
   ! call hdf_write_dataset(group_id, 'v_sl_unscreened', vpuns(:, 1:lmax + 1))
   ! call hdf_write_dataset(group_id, 'v_loc', vp(:, lloc + 1))
   ! call hdf_create_group(group_id, 'v_nl')
   ! call hdf_open_group(group_id, 'v_nl', subgroup_id)
   !
   ! call hdf_close_group(subgroup_id)
   ! call hdf_close_group(group_id)
