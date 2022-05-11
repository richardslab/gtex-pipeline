version 1.0
workflow get_labels_json {
    call do_it
}
task do_it {
    command <<<
        set -euo pipefail

        # Get instance information
        zone=$(curl --silent "http://metadata.google.internal/computeMetadata/v1/instance/zone" -H "Metadata-Flavor: Google")
        instance_name=$(curl --silent "http://metadata.google.internal/computeMetadata/v1/instance/name" -H "Metadata-Flavor: Google")
        
        # Generate an access token.
        # Easier to use other wrappers like $(gcloud auth application-default print-access-token), but this works.
        # Need the SA email used for this VM.
        service_account=$(
            curl \
                --silent \
                --header "Metadata-Flavor: Google" \
                "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/?recursive=true&alt=text" \
            | grep gserviceaccount.com/email \
            | awk '{print $NF}'
        )
        # Can now retrieve an access token.
        access_token=$(
            curl \
                --silent \
                --header "Metadata-Flavor: Google" \
                "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/$service_account/token?alt=text" \
            | grep access_token \
            | awk '{print $NF}'
        )
        
        labels_json=$(
            curl \
            --silent \
            --header "Authorization: Bearer $access_token" \
            https://compute.googleapis.com/compute/v1/${zone}/instances/${instance_name}?fields=labels)

        # Maybe there's a way to get this as something other than json?
        # For now, one will have to parse the results to look for label values.
        echo "$labels_json"
    >>>
    runtime {
        docker: "python:latest"
    }
    output {
        String out = read_string(stdout())
    }
}
