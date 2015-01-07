Vagrant.configure("2") do |config|

    # base box - ubuntu 14.04 lts
    config.vm.box = "ubuntu/trusty64"
    
    # set proper name and ram
    config.vm.provider "virtualbox" do |v|
        v.name = "Algorithmics Vagrantbox"
        v.customize ["modifyvm", :id, "--memory", "2048"]
    end

    # run our provision script to install elasticsearch
    config.vm.provision "shell" do |s|
        s.path = "misc/provision.sh"
    end
end
