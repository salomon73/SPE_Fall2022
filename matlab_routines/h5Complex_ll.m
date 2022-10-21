function z = h5Complex_ll(file, dset)
    fid=H5F.open(file, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
    dset_id=H5D.open(fid, dset);
    dxpl = 'H5P_DEFAULT';
    data = H5D.read(dset_id,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', dxpl);
    z = complex(data.real, data.imaginary);
    H5D.close(dset_id);
    H5F.close(fid);
